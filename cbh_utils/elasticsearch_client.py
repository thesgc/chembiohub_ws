from django.conf import settings
import elasticsearch
import json
import time
try:
    ES_PREFIX = settings.ES_PREFIX
except AttributeError:
    ES_PREFIX = "dev"
ES_MAIN_INDEX_NAME = "chemreg_chemical_index_v2"
from jsonpointer import resolve_pointer

def get_main_index_name():
    return "%s__%s" % (ES_PREFIX, ES_MAIN_INDEX_NAME)


def fix_data_types_for_index( value):
    """Elasticsearch will not index dictionaries"""
    if not value:

        return None
    if isinstance(value, basestring):
        if value.strip():
            return value
        else:
            return None
    if type(value) is dict:
        return json.dumps(value)
    if type(value) is list:
        return [fix_data_types_for_index(v) for v in value]
    return unicode(value)




def build_indexed_fields(document, schema):
    document["indexed_fields"] = []
    for field in schema:
        slashed_json_pointer = "/%s" % field["data"].replace(".", "/")
        raw_value = resolve_pointer(document,slashed_json_pointer, default=None)
        value = fix_data_types_for_index(raw_value)
        if value:
            #We do not add an index for any blank, empty or non existant field, that way
            #we can be sure that the blanks filter will pick up all of the true blank fields
            document["indexed_fields"].append({"name" : field["knownBy"], "value": value, "field_path": field["data"] })

def build_all_indexed_fields(batch_dicts, schema_list):
    assert(len(batch_dicts) == len(schema_list))
    for index,  schema in enumerate(schema_list):
        build_indexed_fields(batch_dicts[index], schema)




def index_dataset(index_name, batch_dicts, schema_list):
    build_all_indexed_fields(batch_dicts, schema_list)
    es_reindex = create_index(
                batch_dicts, index_name)
    if es_reindex.get("errors"):
        print "ERRORS"
        print json.dumps(es_reindex)
        raise Exception("indexing failed")



def create_index(batches, index_name):
    es = elasticsearch.Elasticsearch()
    t = time.time()

    print es.indices.create(
        index_name,
        body=settings.ELASTICSEARCH_INDEX_MAPPING,
        ignore=400)

    bulk_items = []
    if len(batches) > 0:
        for item in batches:
            batch_doc = {
                "update":
                {
                    "_index": index_name,
                    "_type": "newbatches"
                }
            }
            if item.get("id", None):
                batch_doc["update"]["_id"] = str(item["id"])
            bulk_items.append(batch_doc)
            bulk_items.append({"doc" : item, "doc_as_upsert" : True })
        # Data is not refreshed!
        data = es.bulk(body=bulk_items, refresh=True)

    return {}



def get_template_nested_must_clause(field_path, field_query):
    """Match both the original field path and whatever query we are trying to run on that field"""
    column_query = {
                                "term" : {
                                    "indexed_fields.field_path" : field_path
                                }
                            }

    template_must_clause = {
        "nested" : {
                "path" : "indexed_fields",
                "query" : {
                    "bool" : {
                        "must" : [field_query, column_query]
                    }
                }
            }
        }
    return template_must_clause


def build_es_request(queries, textsearch=""):
    must_clauses = []

    if textsearch:
        subquery = {
                    "nested" : {
                            "path" : "indexed_fields",
                            "query" : {
                                
                                    "multi_match" : { 
                                        "type": "phrase_prefix", 
                                        "fields": ["indexed_fields.value",] , 
                                        "query" : textsearch 
                                    }
                                }
                            }
                        }
                    
        must_clauses.append(subquery)
        print subquery


    for query in queries:
        new_query = None
        if query["query_type"] == 'phrase':
            new_query = {
                    "multi_match" :
                    { 
                        "type": "phrase_prefix", 
                        "fields": ["indexed_fields.value",] , 
                        "query" : query["phrase"] 
                    }
                }
        elif query["query_type"] == 'any_of': 
            #We are using a custom char filter but a keyword tokenizer on this field
            #This means that numbers are indexed as their zero-padded equivalents so that they can be sorted in a logical order
            #In order that we can search for a number and have the query also passed through the same filter
             #We need the match query not a terms query
             #The terms query equivalent is shown below. This would work in all cases except 
             #when you are looking for a number.
             #Because the keyword tokenizer does not split up fields into words the match query will always match an exact match just like the terms query does
            # new_query = {
            #         "terms" :
            #         { 
            #             "indexed_fields.value.raw": query["any_of"] 
            #         }
            #     }

            new_query = {                  
                "bool":
                    {
                        "should" : [
                            {"match" : {"indexed_fields.value.raw" : text }} for text in query["any_of"] 
                        ]
                    }
            }

        # elif query["query_type"] == 'starts_with':

        # elif query["query_type"] == 'ends_with':
        elif query["query_type"] ==  'between':
            new_query = {
                    "range" :
                    { 
                        "indexed_fields.value.raw" : 
                            {
                                "gt" : query["greater_than"],
                                "lt"  : query["less_than"]
                            }
                    }
                }
        elif query["query_type"] ==  'greater_than':
            new_query = {
                    "range" :
                    { 
                        "indexed_fields.value.raw" : 
                            {
                                "gt" : query["greater_than"],
                            }
                    }
                }
        elif query["query_type"] ==  'less_than':
            new_query = {
                    "range" :
                    { 
                        "indexed_fields.value.raw" : 
                            {
                                "lt"  : query["less_than"]
                            }
                    }
                }



        if new_query:
            q = get_template_nested_must_clause(query["field_path"], new_query)
            print q
            must_clauses.append(q)
        else:
            new_query =  {
                    "nested" : {
                            "path" : "indexed_fields",
                            "query" : {
                                "term" : {
                                    "indexed_fields.field_path" : query["field_path"]
                                }
                            }
                        }
                    } 
            if query["query_type"] ==  'blanks':
                must_clauses.append({
                    "bool" : {
                            "must_not" : [new_query]
                        }
                    })
            elif query["query_type"] ==  'nonblanks':
                must_clauses.append(new_query)
        

    
            
    return must_clauses


def build_sorts(sorts):
    """This script (written in groovy) picks out 
    the field value in elasticsearch and spits it 
    out as a zero padded string if it is either an integer or a float"""
    elasticsearch_sorts = [
        {
            "_script":{"script":"tmp = ''; for(item in _source.indexed_fields){if(item.field_path==field_path){tmp=item.value}}; if(tmp.isInteger()){tmp = String.format('%014d',tmp.toInteger());};  else if(tmp.isFloat()){def (value1, value2) = tmp.tokenize('.'); tmp = String.format('%014d',value1.toInteger()) + '.' + value2 }; return tmp",
            "params" : {"field_path" : sort["field_path"]},
            "type" : "string", "order" : sort["sort_direction"]}
        }
        for sort in sorts
    ]
    return elasticsearch_sorts



def get_list_data_elasticsearch(queries, index, sorts=[], textsearch="", offset=0, limit=10):
    es = elasticsearch.Elasticsearch()

    if len(queries or textsearch) > 0:
        es_request = {  
                    "query":{
                        
                        "bool" : {
                            "must" : [    
                                   build_es_request(queries, textsearch=textsearch)
                            ]
                        }
                    },
                    "sort" : build_sorts(sorts)
                }
        
    else:
        es_request = {
                        "query" : {"match_all": {}},
                        "sort" : build_sorts(sorts)
                    }
    es_request["from"] = offset
    es_request["size"] = limit

    data = es.search(index, body=es_request)

    return data



def get_detail_data_elasticsearch(index, id):
    es = elasticsearch.Elasticsearch()
    data = es.get(index=index, doc_type="newbatches", id=id)
    return data["_source"]
