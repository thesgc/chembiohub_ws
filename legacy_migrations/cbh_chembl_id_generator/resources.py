from tastypie import fields
from tastypie.resources import Resource, ModelResource
from tastypie.authorization import Authorization
from cbh_chembl_id_generator.models import CBHCompoundId, CBHPlugin
from tastypie.exceptions import ImmediateHttpResponse

def generate_uox_id():
    two_letterg = shortuuid.ShortUUID()
    two_letterg.set_alphabet("ABCDEFGHJKLMNPQRSTUVWXYZ")
    two_letter = two_letterg.random(length=2)
    two_numberg = shortuuid.ShortUUID()
    two_numberg.set_alphabet("0123456789")
    two_number = two_numberg.random(length=2)
    three_letterg = shortuuid.ShortUUID()
    three_letterg.set_alphabet("ABCDEFGHJKLMNPQRSTUVWXYZ")
    three_letter = three_letterg.random(length=3)
    uox_id = "%s%s%s%s" % (
        settings.ID_PREFIX, two_letter, two_number, three_letter)
    try:
        CBHCompoundId.objects.get(assigned_id=uox_id)
        return generate_uox_id()
    except ObjectDoesNotExist:
        return uox_id


class CBHPluginResource(ModelResource):
    handsontable_column = fields.DictField(default=None)

    class Meta:
        always_return_data = True
        queryset = CBHPlugin.objects.all()
        resource_name = 'cbh_plugins'
        authorization = Authorization()
        include_resource_uri = False
        allowed_methods = ['get']
        default_format = 'application/json'

    def dehydrate_handsontable_column(self, bundle):
        return {
            "knownBy": bundle.obj.name,
            "data": "properties.%s" % bundle.obj.space_replaced_name(),
            "readOnly": True,
            "className": "htCenter htMiddle ",
            "renderer": "customFieldRenderer"
        }


class CBHCompoundIdResource(Resource):

    """web service for accessing Compound Ids
    Request without an inchi key will be considered a request for a blinded ID or a forced new ID in that project and installation_key
    Request with and inchi key - inchi will be saved if it is public 
    Inchi will not be saved if it is private

    Scenarios:

    Registering blinded compound in project, need ID
    - asses for uniqueness locally and 
            ask for new ID if new by putting no data in request
            OR ask for new batch (see below)

    New to project Private compound - no data - return new ID
    - Ask for new ID

    Existing in project private or public compound - assigned id- return new batch id
    - Ask for new batch by giving assigned id

    Possibly new Public compound
    - Test for uniqueness against other public compounds - inchi key

    Making private compound public
    - Test for uniqueness against other public compounds - inchi key AND assigned id

    """
    inchi_key = fields.CharField(attribute='inchi_key')
    installation_key = fields.CharField(attribute='installation_key')
    unique_id = fields.CharField(attribute='unique_id')

    class Meta:
        resource_name = 'cbh_compound_ids'
        authorization = Authorization()
        always_return_data = True
        collection_name = "objects"
        allowed_methods = ["post", "patch"]

    def patch_list(self, request, **kwargs):
        """
        Updates a collection in-place.
        The exact behavior of ``PATCH`` to a list resource is still the matter of
        some debate in REST circles, and the ``PATCH`` RFC isn't standard. So the
        behavior this method implements (described below) is something of a
        stab in the dark. It's mostly cribbed from GData, with a smattering
        of ActiveResource-isms and maybe even an original idea or two.
        The ``PATCH`` format is one that's similar to the response returned from
        a ``GET`` on a list resource::
            {
              "objects": [{object}, {object}, ...],
              "deleted_objects": ["URI", "URI", "URI", ...],
            }
        For each object in ``objects``:
            * If the dict does not have a ``resource_uri`` key then the item is
              considered "new" and is handled like a ``POST`` to the resource list.
            * If the dict has a ``resource_uri`` key and the ``resource_uri`` refers
              to an existing resource then the item is a update; it's treated
              like a ``PATCH`` to the corresponding resource detail.
            * If the dict has a ``resource_uri`` but the resource *doesn't* exist,
              then this is considered to be a create-via-``PUT``.
        Each entry in ``deleted_objects`` referes to a resource URI of an existing
        resource to be deleted; each is handled like a ``DELETE`` to the relevent
        resource.
        In any case:
            * If there's a resource URI it *must* refer to a resource of this
              type. It's an error to include a URI of a different resource.
            * ``PATCH`` is all or nothing. If a single sub-operation fails, the
              entire request will fail and all resources will be rolled back.
          * For ``PATCH`` to work, you **must** have ``put`` in your
            :ref:`detail-allowed-methods` setting.
          * To delete objects via ``deleted_objects`` in a ``PATCH`` request you
            **must** have ``delete`` in your :ref:`detail-allowed-methods`
            setting.
        Substitute appropriate names for ``objects`` and
        ``deleted_objects`` if ``Meta.collection_name`` is set to something
        other than ``objects`` (default).
        """
        request = convert_post_to_patch(request)
        deserialized = self.deserialize(request, request.body, format=request.META.get(
            'CONTENT_TYPE', 'application/json'))

        collection_name = self._meta.collection_name

        if collection_name not in deserialized:
            raise BadRequest(
                "Invalid data sent: missing '%s'" % collection_name)

        if len(deserialized[collection_name]) and 'put' not in self._meta.detail_allowed_methods:
            raise ImmediateHttpResponse(response=http.HttpMethodNotAllowed())

        bundles_seen = []

        for data in deserialized[collection_name]:
            # If there's a resource_uri then this is either an
            # update-in-place or a create-via-PUT.

            assigned_id = data.get('assigned_id', '')
            inchi_key = data.get('inchi_key', '')

            if assigned_id and inchi_key:
                data["response"] = CBHCompoundId.objects.make_compound_public(
                    data)
            elif inchi_key:
                data["response"] = CBHCompoundId.objects.new_public_compound(
                    data)
            elif assigned_id:
                data["response"] = CBHCompoundId.objects.new_batch(data)
            else:
                data["response"] = CBHCompoundId.objects.new_private_compound(
                    data)

            bundles_seen.append(data)

        if not self._meta.always_return_data:
            return http.HttpAccepted()
        else:
            to_be_serialized = {}
            to_be_serialized['objects'] = [
                self.full_dehydrate(bundle, for_list=True) for bundle in bundles_seen]
            to_be_serialized = self.alter_list_data_to_serialize(
                request, to_be_serialized)
            return self.create_response(request, to_be_serialized, response_class=http.HttpAccepted)

    def alter_deserialized_list_data(self, request, deserialized):
        schema = deserialized["schema"]
        blinded_key_fieldkeys = [field["key"] for field in schema["form"]]
        for obj in deserialized["objects"]:
            if not obj.get("inchi_key", None) and obj.get("custom_fields"):
                # If no inchi key look for blinded key compounents
                blinded_key_components = [
                    obj["original_installation_keyv"], obj["project_key"]]
                for key in blinded_key_fieldkeys:
                    blinded_key_components += [obj["custom_fields"][key]]
                if len(blinded_key_components) > 2:
                    # If there is a blinded key against the project
                    obj["inchi_key"] = "__".join(blinded_key_components)

        deserialized["project"] = proj
        return deserialized
