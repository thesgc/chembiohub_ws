# -*- coding: utf-8 -*-
import dateutil
import jsonpatch
import requests
import re

'''
A parser object to work with our custom field configs and uncurated fields

'''


class CovertDateOperation(jsonpatch.PatchOperation):

    """Ensures that a data point is formatted correctly for a date field"""

    def apply(self, obj):
        try:
            from_ptr = jsonpatch.JsonPointer(self.location)
        except KeyError as ex:
            raise jsonpatch.InvalidJsonPatch(
                "The operation does not contain a 'from' member")

        subobj, part = from_ptr.to_last(obj)
        try:
            value = subobj[part]
        except (KeyError, IndexError) as ex:
            raise jsonpatch.JsonPatchConflict(str(ex))

        obj = jsonpatch.RemoveOperation({
            'op': 'remove',
            'path': self.location
        }).apply(obj)

        obj = jsonpatch.AddOperation({
            'op': 'add',
            'path': self.location,
            'value': dateutil.parser.parse(value).strftime("%Y-%m-%d")
        }).apply(obj)

        return obj


class SplitOperation(jsonpatch.PatchOperation):

    """Ensures that a data point is formatted correctly for a date field"""

    def apply(self, obj):
        try:
            from_ptr = jsonpatch.JsonPointer(self.location)
        except KeyError as ex:
            raise jsonpatch.InvalidJsonPatch(
                "The operation does not contain a 'from' member")

        subobj, part = from_ptr.to_last(obj)
        try:
            value = subobj[part]
        except (KeyError, IndexError) as ex:
            raise jsonpatch.JsonPatchConflict(str(ex))

        obj = jsonpatch.RemoveOperation({
            'op': 'remove',
            'path': self.location
        }).apply(obj)

        obj = jsonpatch.AddOperation({
            'op': 'add',
            'path': self.location,
            'value': [v.strip() for v in value.split(",") if v.strip()]
        }).apply(obj)

        return obj


class MyJsonPatch(jsonpatch.JsonPatch):

    def __init__(self, patch):
        instance = super(MyJsonPatch, self).__init__(patch)
        self.operations["convertdate"] = CovertDateOperation
        self.operations["split"] = SplitOperation


def apply_json_patch(dictdata, patch):
    mjp = MyJsonPatch(patch)
    try:
        data = mjp.apply(dictdata, in_place=True)
        return data
    except jsonpatch.JsonPatchConflict:
        return dictdata
    except jsonpatch.JsonPointerException:
        return dictdata


def parse_sdf_record(headers, obj, destination_field, mol, fielderrors):
    custom_fields = {}

    for hdr in headers:
        try:
            value = unicode(mol.GetProp(hdr), errors='ignore').strip()
            custom_fields[hdr] = value
            test_specific_parse_errors(hdr, value, fielderrors)

        except KeyError:
            pass
            #custom_fields[hdr]   = u""

    setattr(obj, destination_field, custom_fields)


def parse_pandas_record(headers, obj, destination_field, row, fielderrors, headerswithdata):
    custom_fields = {}

    for hdr in headers:
        if unicode(row[hdr]) == u"nan":
            pass
            # custom_fields[hdr] = ""
        else:
            value = unicode(row[hdr]).strip()
            if value:
                headerswithdata.add(hdr)
                custom_fields[hdr] = unicode(value)
                test_specific_parse_errors(hdr, value, fielderrors)
    # Set excel fields as uncurated
    setattr(obj, destination_field, custom_fields)



def get_uncurated_fields_from_file(correct_file, fielderrors):
    data = correct_file.file.read()
    data = data.replace("\r\n", "\n").replace("\r", "\n")
    ctabs = data.split("$$$$")
    uncurated = []
    for ctab in ctabs:
        pns = re.findall(r'> *<(.+)>',ctab);
        pns2 = re.findall(r'> *<.+> *\S*\n(.+)\n',ctab);
            
        blinded_uncurated_fields = {}
        for idx, hdr in enumerate(pns):
            value = pns2[idx]
            blinded_uncurated_fields[hdr] = unicode(value)
            test_specific_parse_errors(hdr,value, fielderrors)
        uncurated.append(blinded_uncurated_fields)
    return (uncurated, ctabs)

def test_specific_parse_errors(hdr, value, fielderrors):
    if hdr not in fielderrors["stringdate"]:
        try:
            curated_value = dateutil.parser.parse(value).strftime("%Y-%m-%d")
        except:
            fielderrors["stringdate"].add(hdr)
    if hdr not in fielderrors["number"]:
        try:
            curated_value = float(value)
            if hdr not in fielderrors["integer"]:
                if "." in value:
                    fielderrors["integer"].add(hdr)
        except:
            fielderrors["integer"].add(hdr)
            fielderrors["number"].add(hdr)


class APIConverter(object):

    """ writes the formats of APIs to the cbh_core_model format"""

import time


def get_response(uri, session):
    try:
        req = session.get(uri)
        data = req.json()
        for key, value in data.iteritems():
            pass
        return data
    except:
        print data
        print uri
        print "sleeping for 10"
        time.sleep(10)
        return get_response(uri, session)

# https://github.com/akesterson/dpath-python


class ChemblAPIConverter(APIConverter):

    schema_endpoint = "https://www.ebi.ac.uk/chembl/api/data/spore/?format=json"

    def write_schema(self):
        from django.utils.text import capfirst
        from cbh_core_model.models import DataType, PinnedCustomField
        session = requests.Session()
        data = get_response(self.schema_endpoint, session)
        for key, datatype in data["methods"].iteritems():
            if "api_dispatch_list" in key:
                dt, created = DataType.objects.get_or_create(
                    name="chembl__%s" % datatype["collection_name"],
                    uri="https://www.ebi.ac.uk%s/?format=json" % datatype[
                        "schema"],
                    version=data["version"]
                )
                schema = get_response(dt.uri, session)
                for field_name, info in schema.get("fields", {}).iteritems():
                    field_type = info["type"]
                    for fielder, value in PinnedCustomField.FIELD_TYPE_CHOICES.items():
                        if value["data"]["type"] == info["type"]:
                            field_type = fielder
                            break

                    pcf, created = PinnedCustomField.objects.get_or_create(pinned_for_datatype_id=dt.id,
                                                                           name=capfirst(
                                                                               field_name.replace("_", " ")),
                                                                           field_key="/%s/%s/%s" % ("chembl",
                                                                                                    datatype[
                                                                                                        "collection_name"],
                                                                                                    field_name),
                                                                           field_type=field_type,
                                                                           description=info[
                                                                               "help_text"],
                                                                           required=info[
                                                                               "nullable"],
                                                                           position=0)
