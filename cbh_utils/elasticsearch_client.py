"""
This module provides a generic abstraction layer to Elasticsearch
JSON objects can be indexed by providing a schema object per row you are trying to index
Each of the indexed fields will be extracted from the nested JSON using a JSON pointer in the schema object
The actual JSON object will remain unchanged and will be served up by the search functions unchanged.
The fields that you have chosen to index can be queried using JSON Paths in a format specified in the build_es_request function
The web forms required to perform these searches are specified by the JSON schema in the django settings which is compatible
with Angular Schema Form and other web frameworks
This means that an advanced search interface can be provided on top of any JSON-based document data source
Support for foreign key relationships is coming soon
Data from the field will then be prepared for indexing in Elasticsearch

Details of the specific implementation of each step are detailed in the functions 

"""


from django.conf import settings
import elasticsearch
import json
import time
try:
    ES_PREFIX = settings.ES_PREFIX
except AttributeError:
    ES_PREFIX = "dev"
ES_MAIN_INDEX_NAME = "crv3"
from jsonpointer import resolve_pointer
AGG_TERMS_SEPARATOR = "|||"

BATCH_TYPE_NAME = "hashbatches"
import hashlib

ELASTICSEARCH_MAX_FIELD_LENGTH = 256

"""
The mapping patterns were generated using this python function

The reason that we zero pad using regular expressions is to avoid doing that work in python and reduce the size of the request to elasticsearch when indexing

for i in range(1,14):
   data["%ddigit" % i] = {"type": "pattern_replace", "pattern" : "^([1-9]\d{%d})($|\.\d+$)" % (i-1) , "replacement" : "".join(["0" for j in range(0,(13-i))]) + "$1$2"}

The list of patterns to be applied were generated like this:

["%ddigit" % i for i in range(1,14)]
"""


ELASTICSEARCH_INDEX_MAPPING = {
        "settings": {
            "index.store.type": "niofs",
            "index.max_result_window" : 10000000,

            "analysis" : {
                    "char_filter" : {
                        "special_char_space_out" :{ # put spaces around special characters so they can still be indexed
                            "type":"pattern_replace",
                            "pattern":"([()\[\].,\-\+])",
                            "replacement":" $1 "
                        },
                       
                    },
                    "analyzer" : {
                        "default_index" : {
                            "tokenizer" : "whitespace",
                            "filter" : [
                                "lowercase"
                            ],
                            "char_filter" : [
                                "html_strip", "special_char_space_out"
                            ]
                        },
                        "lowercasekeywordanalyzer" : {
                            "tokenizer" : "keyword",
                            "filter" : [
                                "lowercase"
                            ]
                        },
                    }
                },
        },
        "mappings": {
            BATCH_TYPE_NAME: {
                "dynamic": True,
                "_all": {"enabled": False},
                "date_detection": False,
                "properties":{
                    "dataset" : {
                        "type": "object",
                        "dynamic" : False
                    }
                },
                "_source": {
                    "includes": [
                      "*",
                    ],
                    "excludes": [
                      "indexed_fields_*", "sortable_fields_*",
                    ]
                  },
                "dynamic_templates": [
                    {
                    "indexed_fields": {
                        "match": "indexed_fields_*",
                        "match_mapping_type": "string",
                        "mapping": {
                            "type": "string", 
                            "index_options": "positions", 
                            "index": "analyzed", 
                            "omit_norms": True, 
                            "analyzer" : "default_index",
                            "type": "string", 
                            "store": "no", 
                            "fields": {
                                    "raw": {"type": "string", "store": "no", "index": "not_analyzed", "ignore_above": ELASTICSEARCH_MAX_FIELD_LENGTH}
                                }
                            }
                        }
                    },
                    {
                    "sortable_fields": {
                        "match": "sortable_fields_*",
                        "match_mapping_type": "string",
                        "mapping": {
                            "type": "string", "store": "no", "analyzer" : "lowercasekeywordanalyzer", "index": "analyzed", "ignore_above": ELASTICSEARCH_MAX_FIELD_LENGTH
                        }
                    },
                    }

              
               
               
                ]
                
            }
        }
    }

def get_main_index_name():
    return "%s__%s" % (ES_PREFIX, ES_MAIN_INDEX_NAME)

def get_project_index_name(project_id):
    #return get_main_index_name()
    return "%s__%s__project__%d" % (ES_PREFIX, ES_MAIN_INDEX_NAME, project_id)

def get_list_of_indicies(project_ids):
    list_of_indices = [get_project_index_name(project_id) for project_id in project_ids]
    return ",".join(list_of_indices)



def fix_data_types_for_index( value):
    """Elasticsearch will not index dictionaries"""
    if value is None:
        return None

    if not unicode(value):
        #pick up empty strings etc but not false

        return None
    if unicode(value) == "True" or unicode(value) == "False":
        return unicode(value).lower()
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


def get_es_fieldname(pointer):
    m = hashlib.md5()
    m.update(pointer)
    return "indexed_fields_%s" % m.hexdigest()

def get_sortable_es_fieldname(pointer):
    m = hashlib.md5()
    m.update(pointer)
    return "sortable_fields_%s" % m.hexdigest()


def build_indexed_fields(document, schema):
    to_index = {"dataset" : document}
    
    for field in schema:
        slashed_json_pointer = "/%s" % field["data"].replace(".", "/")

        raw_value = resolve_pointer(document,slashed_json_pointer, default=None)
        
        value = fix_data_types_for_index(raw_value)
        if value:
            #We do not add an index for any blank, empty or non existant field, that way
            #we can be sure that the blanks filter will pick up all of the true blank fields
            hashed = get_es_fieldname(field["data"])
            to_index[hashed] = value
            sortable = get_sortable_es_fieldname(field["data"])
            to_index[sortable] = zeropad(value)
    return to_index



def build_all_indexed_fields(batch_dicts, schema_list):
    assert(len(batch_dicts) == len(schema_list))
    for index,  schema in enumerate(schema_list):
        batch_dicts[index] = build_indexed_fields(batch_dicts[index], schema)




def index_dataset( batch_dicts, schema_list, index_names, refresh=True):
    build_all_indexed_fields(batch_dicts, schema_list)
    es_reindex = create_index(
                batch_dicts, index_names, refresh=refresh)
    if es_reindex.get("errors"):
        print "ERRORS"
        print json.dumps(es_reindex)
        raise Exception("indexing failed")



def create_index(batches, index_names, refresh=True):
    es = elasticsearch.Elasticsearch()
    t = time.time()

    for index_name in index_names:
        es.indices.create(
        index_name,
        body=ELASTICSEARCH_INDEX_MAPPING,
        ignore=400)

    bulk_items = []
    result = {}
    if len(batches) > 0:
        for counter, item in enumerate(batches):
            batch_doc = {
                "update":
                {
                    "_index": index_names[counter],
                    "_type": BATCH_TYPE_NAME
                }
            }
            if item["dataset"].get("id", None):
                batch_doc["update"]["_id"] = str(item["dataset"]["id"])
            bulk_items.append(batch_doc)
            bulk_items.append({"doc" : item, "doc_as_upsert" : True })
        # Data is not refreshed!
        result = es.bulk(body=bulk_items, refresh=refresh)

    return result



def build_phase_prefix_query(phrase, field_path):
    return {
                    "multi_match" :
                    { 
                        "type": "phrase_prefix", 
                        "fields": [get_es_fieldname(field_path),] , 
                        "query" : phrase
                    }
                }





def build_es_request(queries, textsearch="", batch_ids_by_project=None):
    must_clauses = []
    if batch_ids_by_project:
        #The postgres backend has converted the chemical search
        #Into a list of ids by project
        #We then join these ids queries on a per project basis
        match_these_ids_by_index = []
        for search_dict in batch_ids_by_project:
            index_names = [get_project_index_name(pid) for pid in search_dict["project_ids"]]
            batch_ids = search_dict["batch_ids"]
            index_query = { "indices":
                                {
                                    "indices": index_names,
                                    "query": {
                                        "ids" : {
                                            "type" : BATCH_TYPE_NAME,
                                            "values" : [str(id) for id in batch_ids]
                                        }
                                    },
                                    "no_match_query" : "none"
                                }
                            }
            match_these_ids_by_index.append(index_query)

        #Each document should be in the specified ID list for the project it is in
        by_index_batch_id_query = {   "bool" :{
                "should" : match_these_ids_by_index
            }
        }
        must_clauses.append(by_index_batch_id_query)



    if textsearch:
        subquery = {"query" : {
                                    "multi_match" : { 
                                        "type": "phrase_prefix", 
                                        "fields": ["indexed_fields_*",] , 
                                        "query" : textsearch 
                                    }
                                }
                            }
                    
        must_clauses.append(subquery)


    for query in queries:
        new_query = None

        if query["query_type"] == 'phrase':
            new_query = build_phase_prefix_query(query["phrase"], query["field_path"])
            
        
        elif query["query_type"] == 'pick_from_list': 

            new_query = {
                    "terms" :
                    { 
                        get_sortable_es_fieldname(query["field_path"] ) : zeropad(query["pick_from_list"])
                    }
                }

            
        elif query["query_type"] ==  'between':
            new_query = {
                    "range" :
                    { 
                         get_sortable_es_fieldname(query["field_path"] ): 
                            {
                                "gt" : zeropad(query["greater_than"]),
                                "lt"  : zeropad(query["less_than"])
                            }
                    }
                }
        elif query["query_type"] ==  'greater_than':
            new_query = {
                    "range" :
                    { 
                         get_sortable_es_fieldname(query["field_path"] ): 
                            {
                                "gt" : zeropad(query["greater_than"]),
                            }
                    }
                }
        elif query["query_type"] ==  'less_than':
            new_query = {
                    "range" :
                    { 
                        get_sortable_es_fieldname(query["field_path"] ) : 
                            {
                                "lt"  : zeropad(query["less_than"])
                            }
                    }
                }
        elif query["query_type"] ==  'nonblanks':
            new_query =  {
                         "exists" : {"field" : "%s.lowercase" % get_es_fieldname(query["field_path"])}
                    } 

        elif query["query_type"] ==  'blanks':
            new_query = {
                    "bool" : {
                            "must_not" : [{
                                     "exists" : {"field" : get_sortable_es_fieldname(query["field_path"] )}
                                } ]
                        }
                    }
        
        if new_query:
            must_clauses.append(new_query)
    
            
    return must_clauses

def get_missing_items_position(direction):
    position_of_missing_items = {
        "asc" : "_last",
        "desc" : "_first"
    }
    return position_of_missing_items[direction]


def build_sorts(sorts):
    """Sort by the lowercase version of the data"""

    elasticsearch_sorts = [
        {
            get_sortable_es_fieldname(sort["field_path"] ) :{
                "order" : sort["sort_direction"],
                "unmapped_type" : "string",
                "missing" : get_missing_items_position(sort["sort_direction"])
            }
        }
        for sort in sorts
    ]
    return elasticsearch_sorts

def get_nested_aggregation_for_field_path(autocomplete_field_path, autocomplete="", autocomplete_size=settings.MAX_AUTOCOMPLETE_SIZE):
    """Based upon an input term and field_path, generate an aggregation to group by that field returning zero padded numbers to get the order right"""
    es_fname = get_es_fieldname(autocomplete_field_path )
    sortable_fname = get_sortable_es_fieldname(autocomplete_field_path )
    base_agg = {
                "field_path_terms" : {
                    "terms" : {
                        "field" : sortable_fname,
                        "size" : autocomplete_size,
                         "order" : { "_term" : "asc" },
                         
                    },
                    "aggs" : {
                                "correct_case":
                                    {
                                        "terms" : {
                                            "field" : "%s.raw" % es_fname,
                                            "size" : 1,
                                             "order" : { "_term" : "asc" },
                                        }
                                    }
                    }
                },
                "unique_count": {
                                    "cardinality" : 
                                    {
                                    "field" : sortable_fname
                                    }
                }
            }

    if autocomplete:
        #If there is a search term, then, having applied a set of filters to the data (project, other search terms)
        #We then try to apply a term filter to the data being aggregated for the search term being looked for
        query =  build_phase_prefix_query(autocomplete, autocomplete_field_path)
    else:
        query = { "match_all": {} }

    base_agg = { 
                "filtered_field_path":
                        {
                        "filter" : query,
                        "aggs" : base_agg

                        }
                }
                    
  
    return base_agg


def remove_existing_queries_for_agg(queries, autocomplete_field_path):
    """We are trying to provide an aggregated autocomplete for a particular field
    Therefore if a query has already been applied to that field, we need to remove that subquery"""
    if autocomplete_field_path:
        new_queries = []
        for q in queries:
            if q["field_path"] != autocomplete_field_path:
                new_queries.append(q)
        return new_queries
    return queries


def get_list_data_elasticsearch(queries, index, sorts=[], autocomplete="", autocomplete_field_path="", autocomplete_size=settings.MAX_AUTOCOMPLETE_SIZE, textsearch="", offset=0, limit=10, batch_ids_by_project=None):
    """Build a query for elasticsearch by going through the input query dictionaries
    Also add autocomplete if required by creating an aggregation"""
    es = elasticsearch.Elasticsearch()
    queries = remove_existing_queries_for_agg(queries, autocomplete_field_path)
    if len(queries) > 0  or len(textsearch) > 0:
        es_request = {  
                    "query":{
                        "indices" : {
                            "indices" : index.split(","), #We ran out of space in a GET request to put all of the indices in, so just using a query indead
                            "query" : {
                                "bool" : {
                                    "filter" : [    
                                           build_es_request(queries, 
                                            textsearch=textsearch, 
                                            batch_ids_by_project=batch_ids_by_project)
                                    ]
                                }
                            },
                            "no_match_query" : "none"
                        }
                    },
                    "sort" : build_sorts(sorts),
                    "_source" : {
                        "include": [ "*" ],
                        "exclude": [ "dataset.bigimage" ]
                    },
                }
        
    else:
        raise Exception("You must input a query")
    if autocomplete_field_path:
        es_request["aggs"] = get_nested_aggregation_for_field_path(autocomplete_field_path, autocomplete=autocomplete, autocomplete_size=autocomplete_size)
    es_request["from"] = offset
    es_request["size"] = limit

    data = es.search("_all", body=es_request, ignore_unavailable=True)

    if autocomplete_field_path:
        for bucket in data["aggregations"]["filtered_field_path"]["field_path_terms"]["buckets"]:
            #Replace the key of the bucket with the first normal case version of it
            for correct_case_bucket in bucket["correct_case"]["buckets"]:
                #If there is a list field there will be choices to choose between - we must pick the
                #Choice where the keys are equal but the correct case bucket may have different case
                if unzeropad(bucket["key"]) == correct_case_bucket["key"].lower():
                    bucket["key"] = correct_case_bucket["key"]

            del bucket["correct_case"]

    return data

def zeropad(input_string):
    """Zero pad integers and floats to ensure that they can be sorted correctly"""
    if isinstance(input_string, basestring):
        if len(input_string) > 0 :
            if input_string[0].isdigit():
                if "." in input_string:
                    if input_string.replace(".", "", 1).isdigit():
                        data = input_string.split(".")
                        
                        return "%s.%s" % (data[0].zfill(13), data[1])
                elif input_string.isdigit():
                    return input_string.zfill(13)
        return input_string.lower()
    elif hasattr(input_string, '__iter__') :
        """we expect a list here"""
        return [zeropad(string) for string in input_string]
    else:
        return input_string
    



def unzeropad(input_string):
    """Remove zero padding from integers and floats"""
    replace_up_to = 0
    if len(input_string) > 0:
        if input_string[0] == "0":
            if input_string.replace(".", "", 1).isdigit():
                for index, char in enumerate(input_string):
                    if char != "0":
                        replace_up_to = index
                        #found the first non zero so break
                        break 
    return input_string[replace_up_to:]


def get_detail_data_elasticsearch(index, id):
    es = elasticsearch.Elasticsearch()
    es_request = {
                        "query":{
                            "indices" : {
                                "indices" : index.split(","), #We ran out of space in a GET request to put all of the indices in, so just using a query indead
                                "query" : {
                                    "term": {"_id" : id},
                                },
                                "no_match_query" : "none"
                            }
                        },


                        "_source" : {
                            "include": [ "*" ],
                            "exclude": [ "indexed_fields.*" ]
                        }
                    }
    data = es.search("_all", body=es_request, ignore_unavailable=True)

    return data["hits"]["hits"][0]["_source"]["dataset"]
