# -*- coding: utf-8 -*-
import dateutil
import jsonpatch
import requests
import re

'''
A parser object to work with our custom field configs and uncurated fields

'''

def get_all_sdf_headers(filename):
    """Use the unix tools grep, cut, sort and uniq to pull out a set of headers from a file"""
    from subprocess import Popen, PIPE
    from shlex import split
    p1 = Popen(split('grep "^>" "%s"' % filename), stdout=PIPE)
    p2 = Popen(split('cut -d "<" -f2'), stdin=p1.stdout, stdout=PIPE)
    p3 = Popen(split('cut -d ">" -f1'), stdin=p2.stdout, stdout=PIPE)
    p4 = Popen(split('sort'), stdin=p3.stdout, stdout=PIPE)
    p5 = Popen(split('uniq'), stdin=p4.stdout, stdout=PIPE)
    out = p5.communicate()
    return [i for i in out[0].split("\n") if i]
    
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

def get_sdf_count(correct_file):

    from subprocess import Popen, PIPE
    from shlex import split
    p1 = Popen(split('grep -c "\$\$\$\$" "%s"' % correct_file.file.name), stdout=PIPE)
    out = p1.communicate()
    print out
    return int(out[0].strip())
    

def get_uncurated_fields_from_file(correct_file, fielderrors):
    data = correct_file.file.read().decode('string-escape').decode("utf-8", "ignore")
    data = data.replace("\r\n", "\n").replace("\r", "\n")
    ctabs = data.split("$$$$\n")
    uncurated = []
    length = len(ctabs)
    ctab_parts = []
    for index, ctab in enumerate(ctabs):
        pns = re.findall(r'> *<(.+)> *\S*\n+(.+)\n+',ctab)
        blinded_uncurated_fields = {}
        for key, value in pns:
            blinded_uncurated_fields[key] = value
            test_specific_parse_errors(key,value, fielderrors)

        uncurated.append(blinded_uncurated_fields)
        if index  == length -1:
            
            if ctab.strip():
                # The last item may not have a line break
                # And therefore the dollars must be removed
                ctabs[index] = ctab.split("$$$$")[0]
            else:
                # If the last item is empty it is not 
                # a real ctab so skip it
                ctabs = ctabs[:-1]
                continue

        ctab_part = ctab.split("END\n")[0] + "END"
        ctab_parts.append(ctab_part)

    return (uncurated, ctabs, ctab_parts)

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