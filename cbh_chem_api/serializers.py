# -*- coding: utf-8 -*-
"""
Original serializer module for CBHCompoundbatchResource
Currently in the process of being deprecated"""
import cStringIO

from tastypie.serializers import Serializer

import re
import json
import pandas as pd


import base64
import copy
import os

from cbh_core_api.resources import get_field_name_from_key
from cbh_core_api.resources import get_key_from_field_name
from cbh_core_api import parser as coreparser

from tastypie.exceptions import ImmediateHttpResponse, BadRequest

def flatten_dict(d, base=None):
    """Converts a dictionary of dictionaries or lists into a simple
    dictionary.

    For example the following dictionary

    foobar = {'key1': 'value1',
              'key2': {'skey1': 'svalue1'},
              'key3': ['svalue2', 'svalue3']}

    gets converted to

    foobar = {'key1': 'value1',
              'key2.skey1': 'svalue1',
              'key3.0': 'svalue2',
              'key3.1': 'svalue3'}

    """
    new_dict = {}
    for key, value in d.iteritems():
        if isinstance(value, dict):
            new_base = ''
            if base:
                new_base = '%s.' % base
            new_base += key
            new_dict.update(flatten_dict(value, base=new_base))
        elif isinstance(value, list):
            new_base = ''
            if base:
                new_base += '%s.' % base
            new_base += '%s' % key
            i = 0
            for item in value:
                new_base_index = new_base + '.%d' % i
                if isinstance(item, dict):
                    new_dict.update(flatten_dict(item, base=new_base_index))
                else:
                    new_dict.update({new_base_index: item})
                i += 1
        elif base:
            new_dict.update({'%s.%s' % (base, key): value})
        else:
            new_dict.update({key: value})

    return new_dict



class CBHCompoundBatchSerializer(Serializer):
    # formats = ['json']
    # content_types = {
    #     'json': 'application/json',
    # }

    def to_json(self, data, options=None):
        # Changes underscore_separated names to camelCase names to go from
        # python convention to javacsript convention
        data = self.to_simple(data, options)

        def underscoreToCamel(match):
            return match.group()[0] + match.group()[2].upper()

        def camelize(data):
            if isinstance(data, dict):
                new_dict = {}
                for key, value in data.items():
                    new_key = re.sub(r"[a-z]_[a-z]", underscoreToCamel, key)
                    if new_key in ["customFields", "uncuratedFields"]:

                        for k, v in value.iteritems():
                            if isinstance(v, basestring):
                                if v.startswith("[") and v.endswith("]"):
                                    try:
                                        value[k] = json.loads(v)
                                        continue
                                    except:
                                        pass
                                elif "." in v:
                                    try:
                                        value[k] = float(v)
                                        continue
                                    except:
                                        pass
                                else:
                                    try:
                                        value[k] = int(v)
                                        continue
                                    except:
                                        pass
                                value[k] = v
                        new_dict[new_key] = value
                    elif new_key in ["projectfull"]:
                        new_dict[key] = value

                    else:
                        new_dict[new_key] = camelize(value)
                return new_dict
            if isinstance(data, (list, tuple)):
                for i in range(len(data)):
                    data[i] = camelize(data[i])
                return data
            return data

        camelized_data = camelize(data)
        for key, value in camelized_data.iteritems():
            try:
                dictd = json.loads(value)
                if isinstance(dictd, dict):
                    camelized_data[key] = dictd
            except:
                pass
        return json.dumps(camelized_data, sort_keys=True)

    def from_json(self, content):
        # Changes camelCase names to underscore_separated names to go from
        # javascript convention to python convention
        data = json.loads(content)

        def camelToUnderscore(match):
            return match.group()[0] + "_" + match.group()[1].lower()

        def underscorize(data):
            if isinstance(data, dict):
                new_dict = {}
                for key, value in data.items():
                    new_key = re.sub(r"[a-z][A-Z]", camelToUnderscore, key)
                    if new_key in ["custom_fields", "uncurated_fields", "compoundstats", "batchstats", "warnings", "properties"]:
                        new_dict[new_key] = value

                    else:
                        new_dict[new_key] = underscorize(value)
                return new_dict
            if isinstance(data, (list, tuple)):
                for i in range(len(data)):
                    data[i] = underscorize(data[i])
                return data
            return data

        underscored_data = underscorize(data)

        return underscored_data




def convert_query(data):
    if isinstance(data, dict):
        new_dict = {}
        for key, value in data.items():
            new_key = get_key_from_field_name(key)
            new_key = new_key.replace("uncuratedFields", "uncurated_fields")
            new_key = new_key.replace("customFields", "custom_fields")

            new_dict[new_key] = convert_query(value)
        return new_dict
    if isinstance(data, (list, tuple)):
        for i in range(len(data)):
            data[i] = convert_query(data[i])
        return data
    return data


def whitespaced(string):
    if string:
        s = re.sub('[^0-9a-zA-Z]+', ' ', unicode(string))
        return ' %s' % s
    else:
        return ""


def get_agg(field_name, field_value):
    return '%s|%s|%s|%s' % (field_name,
                            field_value,
                            whitespaced(field_name),
                            whitespaced(field_value))


class CBHCompoundBatchElasticSearchSerializer(Serializer):
    formats = ['json']
    content_types = {
        'json': 'application/json',
    }

    def convert_query(self, es_request):

        def camelToUnderscore(match):
            return match.group()[0] + "_" + match.group()[1].lower()

        es_request["query"] = convert_query(es_request["query"])
        es_request["sort"] = convert_query(es_request["sort"])
        newsort = []
        for item in es_request["sort"]:
            newItem = {}
            for sort, direction in item.items():
                if ("." in sort):
                    newItem[sort + "___sortable"] = direction
                else:
                    new_key = re.sub(r"[a-z][A-Z]", camelToUnderscore, sort)
                    if new_key != "id":
                        new_key += ".raw"
                    newItem[new_key] = direction
            newsort.append(newItem)
        es_request["sort"] = newsort

    def handle_data_from_django_hstore(self, value):
        '''Hstore passes data in the wrong format'''
        for k, v in value.iteritems():
            if isinstance(v, basestring):
                if v.startswith("[") and v.endswith("]"):
                    try:
                        value[k] = json.loads(v)
                        continue
                    except:
                        pass
                # elif "." in v:
                #     try:
                #         value[k] = float(v)
                #         continue
                #     except:
                #         pass
                # else:
                #     try:
                #         value[k] = int(v)
                #         continuecompound_stats
                #     except:
                #         pass
                value[k] = unicode(v)




    def to_es_ready_data(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        data['custom_field_list'] = []
        self.handle_data_from_django_hstore(data["custom_fields"])

        for key, value in data["custom_fields"].items():
            if not key.endswith("___sortable"):
                if type(value) == list:
                    for val in value:
                        if val:
                            val = val.replace(u"\n|\r", " ")
                            data['custom_field_list'].append(
                                {'name': key,
                                 'value': val,
                                 'searchable_name': key.split(" ")[0].lower(),
                                 'aggregation': get_agg(key, val)}
                            )
                else:
                    if value:
                        value = value.replace(u"\n|\r", " ")
                        data['custom_field_list'].append(
                            {'name': key,
                             'value': value,
                             'searchable_name': key.split(" ")[0].lower(),
                             'aggregation': get_agg(key, value)})

        for key, value in data.items():
            if key in ["custom_fields", "uncurated_fields"]:
                if options and options.get("underscorize", False):
                    data[key] = self.underscorize_fields(value)
                self.handle_data_from_django_hstore(value)
        return data

    def to_es_ready_non_chemical_data(self, data, options=None):
        options = options or {}
        newdata = {}

        data = self.to_simple(data, options)

        newdata = copy.deepcopy(data)
        newdata['custom_field_list'] = []
        self.handle_data_from_django_hstore(data["custom_fields"])

        for key, value in data["custom_fields"].items():
            if not key.endswith("___sortable"):
                if type(value) == list:
                    for val in value:
                        if val:
                            v = val.replace('\n', ' ').replace('\r', '')
                            agg = '%s|%s' % (key, v)
                            newdata['custom_field_list'].append(
                                {'name': key, 'value': v, 'searchable_name': key.split(" ")[0].lower(), 'aggregation': agg})
                else:
                    if value:
                        v = value.replace('\n', ' ').replace('\r', '')
                        agg = '%s|%s' % (key, v)

                        newdata['custom_field_list'].append(
                            {'name': key, 'value': v, 'searchable_name': key.split(" ")[0].lower(), 'aggregation': agg})

        # newdata['custom_field_list'].append({'name': "Project", 'value':newdata['project'], 'searchable_name': 'project', 'aggregation': '%s|%s' % ('Project', newdata['project']) })
        # newdata['custom_field_list'].append({'name': "Upload Id", 'value':newdata['multiple_batch_id'], 'searchable_name': 'upload', 'aggregation': '%s|%d' % ('Upload', newdata['multiple_batch_id']) })

        for key, value in data.items():
            if key in ["custom_fields", "uncurated_fields"]:
                if options and options.get("underscorize", False):
                    newdata[key] = self.underscorize_fields(value)
                self.handle_data_from_django_hstore(value)
        return newdata

    def to_python_ready_data(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        data.pop("custom_field_list", False)
        for key, value in data.items():
            if key in ["custom_fields", "uncurated_fields"]:
                data[key] = self.deunderscorize_fields(value)
        return data

    def to_json(self, data, options=None):
        self.to_es_ready_data(data, options=options)
        return json.dumps(data, sort_keys=True)

    def make_sortable_data(self, value):
        "Zero pad an integers or floats"
        if type(value) is list:
            return [self.make_sortable_data(item) for item in value]
        if isinstance(value, basestring):
            value = value.strip()
            if value.replace(".", "", 1).isdigit():
                splitup = value.split(".")
                if len(splitup) == 2:
                    return "%s.%s" % (splitup[0].zfill(14) , splitup[1] )
                return value.zfill(14)
        return value



    def underscorize_fields(self, dictionary):
        data = {
            get_key_from_field_name(key): value
            for key, value in dictionary.items()
        }
        for key, value in data.items():
            data[key + "___sortable"] = self.make_sortable_data(value)
        return data


    def deunderscorize_fields(self, dictionary):
        return {
            get_field_name_from_key(key): value
            for key, value in dictionary.items() if not key.endswith("___sortable")
        }


