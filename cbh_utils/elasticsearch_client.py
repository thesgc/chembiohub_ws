from django.conf import settings
import elasticsearch
import json
import time
try:
    ES_PREFIX = settings.ES_PREFIX
except AttributeError:
    ES_PREFIX = "dev"
ES_MAIN_INDEX_NAME = "chemreg_chemical_index_v2"


def get_main_index_name():
    return "%s__%s" % (ES_PREFIX, ES_MAIN_INDEX_NAME)



def create_temporary_index(batches, request, index_name):
    es = elasticsearch.Elasticsearch()
    t = time.time()
    store_type = "niofs"

    create_body = {
        "settings": {
            "index.store.type": store_type,
            "analysis" : {
                    "analyzer" : {
                        "default_index" : {
                            "tokenizer" : "whitespace",
                            "filter" : [
                                "lowercase"
                            ],
                            "char_filter" : [
                                "html_strip"
                            ]
                        }
                    }
                },
        },
        "mappings": {
            "_default_": {
                "_all": {"enabled": False},
                "date_detection": False,
                
                "properties":{
                    "indexed_fields" :{
                        "type": "nested",
                        "properties" : {
                                "name": {
                                    "type": "string", 
                                    "index": "not_analyzed"
                                },
                                "value":  
                                      {
                                        "type": "string", 
                                        "store": "no", 
                                        "index_options": "positions", 
                                        "index": "analyzed", 
                                        "omit_norms": True, 
                                        "analyzer" : "default_index",
                                        "fields": {
                                            "raw": {"type": "string", "store": "no", "index": "not_analyzed", "ignore_above": 256}
                                        }
                                    }
                            }
                        }
                },
             "dynamic_templates": [{
                    "ignored_fields": {
                        "match": "*",
                        "match_mapping_type": "string",
                        "mapping": {
                            "type": "string", "store": "no", "include_in_all": False, "index" : "no"
                        }
                    }
                }]
                
            }
        }
    }

    es.indices.create(
        index_name,
        body=create_body,
        ignore=400)

    bulk_items = []
    if len(batches) > 0:
        for item in batches:
            bulk_items.append({
                "index":
                {
                    "_id": str(item["id"]),
                    "_index": index_name,
                    "_type": "batches"
                }
            })
            bulk_items.append(item)
        # Data is not refreshed!
        return es.bulk(body=bulk_items, refresh=True)
    return {}
