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
from jsonschema import validate

def get_main_index_name():
    return "%s__%s" % (ES_PREFIX, ES_MAIN_INDEX_NAME)


def fix_data_types_for_index( value):
    """Elasticsearch will not index dictionaries"""
    if not value:
        return "__EMPTY"
    if isinstance(value, basestring):
        return value
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
                "index":
                {
                    "_index": index_name,
                    "_type": "newbatches"
                }
            }
            if item.get("id", None):
                batch_doc["_id"] = str(item["id"])
            bulk_items.append(batch_doc)
            bulk_items.append(item)
        # Data is not refreshed!
        data = es.bulk(body=bulk_items, refresh=True)
        print data
        return data
    return {}

def get_schema_for_query_type(query_type):
    if not query_type:
        raise Exception("Query Type must not be null, cannot continue")
    schema = None
    for qt in settings.CBH_QUERY_TYPES:
        if qt["value"] == query_type:
            required_fields = settings.CBH_QUERY_ALWAYS_REQUIRED + qt["required_fields"]
            print required_fields
            schema = { req :settings.CBH_QUERY_SCHEMA[req] for req in required_fields}
    if not query_type:
        raise Exception("Query Type must be in available query types, cannot continue")
    return schema

def validate_elasticsearch_queries(queries):
    for query in queries:
        qt = query.get("query_type", None)
        schema = get_schema_for_query_type(qt)
        validate(query, schema)


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


def build_es_request(queries):
    must_clauses = []
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

    base_query = {  
                    "query":{
                        
                        "bool" : {
                            "must" : [    
                                   must_clauses
                            ]
                        }
                    }
                }
            
        
    return base_query


def get_list_data_elasticsearch(queries, index, offset=0, limit=10):
    es = elasticsearch.Elasticsearch()
    validate_elasticsearch_queries(queries)
    if len(queries) > 0:
        es_request = build_es_request(queries)
    else:
        es_request = {"query" : {"match_all": {}}}
    es_request["from"] = offset
    es_request["size"] = limit

    data = es.search(index, body=es_request)

    return data
