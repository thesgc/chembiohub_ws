"""
deprecated soon - the original elasticsearch client for compound batches, now only used for the indexing of
new data in a temporary index
"""
from django.conf import settings
import elasticsearch
import json
import time
try:
    ES_PREFIX = settings.ES_PREFIX
except AttributeError:
    ES_PREFIX = "dev"
ES_MAIN_INDEX_NAME = "chemreg_chemical_index"


def get_temp_index_name(request, multi_batch_id):
    """Generate an index name for a particular multiple batch"""
    index_name = "%s__temp_multi_batch__%s__%s" % (
        ES_PREFIX, request.session.session_key, str(multi_batch_id))
    return index_name
    # return "%s__temp_multi_batch__%s__%s" % (ES_PREFIX,
    # request.session.session_key, str(multi_batch_id))


def get_main_index_name():
    """Generate the main index name for elasticsearch"""
    return "%s__%s" % (ES_PREFIX, ES_MAIN_INDEX_NAME)


def delete_index(index_name):
    """Delete an index permanently by name"""
    es = elasticsearch.Elasticsearch()
    return es.indices.delete(index_name,  ignore=[400, 404])


def get_action_totals(index_name,  bundledata):
    """
    Use elasticsearch aggregations to find the stats required for the preview 
    UI when uploading a set of compounds or inventory items
    """
    es_request_body = {
        "size": 0,
        "query": {"match_all": {}},
        "aggs": {
            "actions": {
                "terms": {"field": "properties.action.raw"}
            }
        }
    }
    es = elasticsearch.Elasticsearch()

    result = es.search(index_name, body=es_request_body)
    bundledata["savestats"] = {"ignoring": 0, "newbatches": 0}
    for buck in result["aggregations"]["actions"]["buckets"]:
        if buck.get("key", "") == "New Batch":
            bundledata["savestats"]["newbatches"] = buck.get("doc_count", 0)
        if buck.get("key", "") == "Ignore":
            bundledata["savestats"]["ignoring"] = buck.get("doc_count", 0)
    return bundledata


def get(index_name, es_request_body, bundledata):
    """Perform a search request against the elasticsearch index specified"""
    es = elasticsearch.Elasticsearch()
    result = es.search(index_name, body=es_request_body)
    data = []
    for hit in result["hits"]["hits"]:
        data.append(hit["_source"])

    bundledata["meta"] = {"totalCount": result["hits"]["total"]}
    bundledata["objects"] = data
    if result.get("aggregations", None):
        bundledata["aggregations"] = [res["key"]
            for res in result["aggregations"]["autocomplete"]["buckets"]]

    return bundledata

def get_project_uri_terms(project_id_list):
    """Because of the way that we indexed the data in ChemiReg V1 the project is represented as a URI
    but the actual search requests and the permissions requests tended to contain list oif ids
    This function converts a list of ids into a list of resource uris"""
    project_terms = []
    for proj in project_id_list:
        project_name = '/%s/cbh_projects/%d' % (
            settings.WEBSERVICES_NAME, proj)
        project_terms.append({'term': {'project.raw': project_name}})
    return [{'bool': {
        'should': project_terms,
    }, }]



def get_custom_fields_query_from_string(cf_string):
    """deperecated"""
    project_terms = []
    for keyvalue in json.loads(cf_string):
        
        splitted = keyvalue.split("|")
        project_terms.append({"term": 
                {"customFields.%s.raw" % splitted[0] : splitted[1]}}
                )
    return {'bool': {
            'should': project_terms,
                }, }



def get_cf_aggregation(search_term, field, single_field):
    """Deperecated way of building an aggregation for a custom field"""
    search_regex = '.*%s.*|.*%s.*|.*%s.*|.*%s.*' % (
        search_term.title(), search_term, search_term.upper(), search_term.lower())
    field_to_search = '%s.raw' % (field)
    should_filter = [
                         
                         {'regexp': {
                             'custom_field_list.value.raw':  search_regex}}
                         ]
    if single_field:
        #We limit the search to a single field
        search_regex = '^%s.*' % (single_field, )
    else:
        #We search only 
        should_filter += [{'prefix': {
                             'custom_field_list.searchable_name.raw':  search_term.lower()}}]


    search_filter =  {'bool': {
            'should': should_filter
                    }, 
                    }

    aggs = {
            'autocomplete': {
                'terms': {'field': field_to_search,
                          'size': 300,
                          'include': str(search_regex)
                          }
            }
        }
    return (search_filter, aggs)



def get_autocomplete(projects, search_term, field, custom_fields=None, single_field=None):
    """deprecated autocomplete method"""
    project_terms = get_project_uri_terms(projects)
    es = elasticsearch.Elasticsearch()
    
    # Search for a space before item to ensure it is a separate word, note
    # that the terms have been formatted to allow this type of match to work
    search_regex = '.*%s.*|.*%s.*|.*%s.*|.*%s.*' % (
        search_term.title(), search_term, search_term.upper(), search_term.lower())
    field_to_search = '%s.raw' % (field)
   

    must_list = project_terms

    if search_term and custom_fields:
        must_list.append({'bool': {
            'should': [
                         {'prefix': {
                             'custom_field_list.searchable_name.raw':  search_term.lower()}},
                         {'regexp': {
                             'custom_field_list.value.raw':  search_regex}}
                         ],
        }, })


    if (custom_fields and single_field):
        # create a bool must term which is the custom field identifier
        #cust_str = 'custom_fields.value.raw' % (single_field)
        # agg_regex = '^%s(.*)(%s)' % (single_field, search_regex)
        agg_regex = '^%s.*' % (single_field, )

    else:
        agg_regex = search_regex
        # must_list.append({
        #                       'bool': {
        #                           'must': {'prefix': { 'custom_field_list.name.raw':  single_field } }
        #                       }
        #                     })
    must_list.append({
        'filtered': {
                'filter' : {"bool":
                                   {"should": [{"term": {"properties.archived": "false"}},
                                               {"missing": {"field": "properties.archived"}}]}
                                   },
            }
        } )

    body = {
        'query': {
            'bool': {
                'must': must_list
            },
        },
        'aggs': {
            'autocomplete': {
                'terms': {'field': field_to_search,
                          'size': 300,
                          'include': str(agg_regex)
                          }
            }
        },
        
        'size': 0,
    }

    result = es.search(body=body)
    # return the results in the right format
    data = [res["key"]
            for res in result["aggregations"]["autocomplete"]["buckets"]]
    return data


def create_temporary_index(batches, request, index_name):
    """Index data in the old format, now mostly deprecated except for the file upload preview"""
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
                    "created": {
                        "type" : "string",
                        "index": "not_analyzed"
                    },
                },
                "dynamic_templates": [{
                    "ignored_fields": {
                        "match": "ctab|std_ctab|canonical_smiles|original_smiles",
                        "match_mapping_type": "string",
                        "mapping": {
                            "type": "string", "store": "no", "include_in_all": False, "index" : "no"
                        }
                    }
                },
                {
                    "sortable": {
                        "match": "customFields.*___sortable",
                        "match_mapping_type": "string",
                        "mapping":  {
                            "type": "string", 
                            "store": "no", 
                            "index": "not_analyzed", 
                            "ignore_above": 256
                        }
                    }
                },
                {
                    "uncurated_fields": {
                        "match": "uncuratedFields.*|image|bigimage",
                        "match_mapping_type": "string",
                        "mapping": {
                            "type": "string", "index": "no", "include_in_all": False
                        }
                    }
                },
                {
                    "string_fields": {
                        "match": "*",
                        "match_mapping_type": "string",
                        "mapping": {
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
                ]
            }
        }
    }
    # if(index_name == get_main_index_name()):
    #     create_body['mappings']['_source'] = { 'enabled':False }
    #index_name = get_temp_index_name(request, multi_batch_id)

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


def get_project_index_name(project):
    index_name = "%s__project__%s" % (ES_PREFIX, str(project.id))
    return index_name


def reindex_compound(dataset, id):
    """Reindex an item now deprecated"""
    # reindex the specified compound in the specified index
    index_name = get_main_index_name()
    es = elasticsearch.Elasticsearch()
    update_body = dataset
    return es.index(id=id, doc_type="batches", index=index_name, body=update_body, refresh=True)
