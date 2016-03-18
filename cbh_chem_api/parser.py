


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
