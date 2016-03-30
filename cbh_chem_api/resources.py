from chembl_business_model.models.compounds import CompoundProperties, MoleculeDictionary
from cbh_chembl_model_extension.models import CBHCompoundBatch, generate_uox_id
from tastypie.resources import ALL
from tastypie.resources import ALL_WITH_RELATIONS
from tastypie.resources import ModelResource, Resource
from tastypie.serializers import Serializer
from django.conf import settings
from django.conf.urls import *
from django.core.exceptions import ObjectDoesNotExist
from tastypie.authorization import Authorization
from tastypie import fields
from cbh_core_api.authorization import ProjectAuthorization
from cbh_core_api.resources import SimpleResourceURIField, UserResource, UserHydrate, CBHNoSerializedDictField,  ChemregProjectResource, ChemRegCustomFieldConfigResource, NoCustomFieldsChemregProjectResource
from cbh_utils import elasticsearch_client
import json 
from django.http import HttpResponse, HttpRequest
from tastypie import http
from tastypie.exceptions import Unauthorized, BadRequest
from tastypie.utils.mime import determine_format, build_content_type
from base64 import b64decode, b64encode
import time
from cbh_chem_api.new_serializers import CBHCompoundBatchSerializer
from django_q.tasks import async_iter, result

EMPTY_ARRAY_B64 = b64encode("[]")

class CompoundPropertiesResource(ModelResource):
    class Meta:
        queryset = CompoundProperties.objects.all()
        fields = ["alogp", "full_molformula", "full_mwt"]

class MoleculeDictionaryResource(ModelResource):
    compoundproperties = fields.ForeignKey(CompoundPropertiesResource, 'compoundproperties',  null=True, readonly=True, full=True)
    authorization = ProjectAuthorization()
    project = fields.ForeignKey(
        ChemregProjectResource, 'project', blank=False, null=False)


    class Meta:
        queryset = MoleculeDictionary.objects.all()
        fields = ["compoundproperties"]



def reformat_project_data_fields_as_table_schema( table_schema_type, project_data_fields_json):
    """Takes the standard output from project data fields and reformats it for use in table formats 
    by including static schemas from the application settings - these include the fields that should be listed before and after the given field"""
    if table_schema_type not in settings.TABULAR_DATA_SETTINGS:
        raise Exception("Trying to reformat for a format that is not configured in TABULAR_DATA_SETTINGS")
    start_items = settings.TABULAR_DATA_SETTINGS[table_schema_type]["start"]
    end_items = settings.TABULAR_DATA_SETTINGS[table_schema_type]["end"]
    start_schema = [settings.TABULAR_DATA_SETTINGS["schema"][start] for start in start_items]
    middle_table_schema = [field["edit_form"]["form"][0] for field in project_data_fields_json]
    end_schema = [settings.TABULAR_DATA_SETTINGS["schema"][end] for end in end_items]

    return start_schema + middle_table_schema + end_schema


class CBHChemicalSearch(Resource):
    uuid = fields.CharField()
    query_type = fields.CharField()
    molfile = fields.CharField()
    smiles = fields.CharField()
    pids = fields.ListField()
    task_id = fields.CharField()



    def post_list(self, request, **kwargs):
        """
        Creates a new resource/object with the provided data.
        Calls ``obj_create`` with the provided data and returns a response
        with the new resource's location.
        If a new resource is created, return ``HttpCreated`` (201 Created).
        If ``Meta.always_return_data = True``, there will be a populated body
        of serialized data.
        """
        deserialized = self.deserialize(request, request.body, format=request.META.get('CONTENT_TYPE', 'application/json'))
        deserialized = self.alter_deserialized_detail_data(request, deserialized)
        bundle = self.build_bundle(data=dict_strip_unicode_keys(deserialized), request=request)
        updated_bundle = self.obj_create(bundle, **self.remove_api_resource_names(kwargs))
        
        updated_bundle = self.alter_detail_data_to_serialize(request, updated_bundle)
        return self.create_response(request, updated_bundle, response_class=http.HttpCreated, location=None)


    def obj_create(self, bundle, **kwargs):
        """
        Generate a task id and send the relevant tasks to django q
        """
        pids = bundle.data["pids"]
        args = [(int(pid), bundle.data["smiles"]) for pid in pids.split(",")]
        bundle.data = self.add_task_id(bundle.data, args)
        return bundle
        
    def add_task_id(self, data, args):
        data["task_id"] = async_iter('cbh_chem_api.tasks.get_structure_search_for_project', args)
        return data




class BaseCBHCompoundBatchResource(UserHydrate, ModelResource):
    uuid = fields.CharField(default="")
    timestamp = fields.CharField(default="")
    project = fields.ForeignKey(
        ChemregProjectResource, 'project', blank=False, null=False)
    creator = SimpleResourceURIField(UserResource, 'created_by_id', null=True, readonly=True)
    projectfull = fields.ForeignKey(
         ChemregProjectResource, 'project', blank=False, null=False, full=True, readonly=True)
    related_molregno = fields.ForeignKey(MoleculeDictionaryResource, 'related_molregno',  null=True, readonly=True, full=True)
    uncurated_fields = CBHNoSerializedDictField('uncurated_fields')
    warnings = CBHNoSerializedDictField('warnings')
    properties = CBHNoSerializedDictField('properties')
    custom_fields = CBHNoSerializedDictField('custom_fields')



    def dehydrate_timestamp(self, bundle):
        return str(bundle.obj.created)[0:10]

    def dehydrate_uuid(self, bundle):
        """This is either a blinded batch or a uuid"""
        
        if bundle.obj.related_molregno:
            if bundle.obj.related_molregno.chembl:
                if bundle.obj.related_molregno.chembl.chembl_id:
                    return bundle.obj.related_molregno.chembl.chembl_id
        if bundle.obj.blinded_batch_id.strip():
            return bundle.obj.blinded_batch_id

    def create_response(self, request, data, response_class=HttpResponse, **response_kwargs):
        """
        Extracts the common "which-format/serialize/return-response" cycle.
        Mostly a useful shortcut/hook.
        """
        desired_format = self.determine_format(request)
        if response_class == http.HttpCreated or response_class == http.HttpAccepted:
            batches = []
            if data.data.get("objects"):
                for b in data.data.get("objects"):
                    batches.append(b.obj)
            else:
                batches.append(data.obj)
            index_batches_in_new_index(batches)

        serialized = self.serialize(request, data, desired_format)
        rc = response_class(content=serialized, content_type=build_content_type(
            desired_format), **response_kwargs)
        if(desired_format == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'):
            rc['Content-Disposition'] = 'attachment; filename=export_from_chembiohub_chemreg%d.xlsx' % int(time.time())
        elif(desired_format == 'chemical/x-mdl-sdfile'):
            rc['Content-Disposition'] = 'attachment; filename=export_from_chembiohub_chemreg%d.sdf' % int(time.time())
        return rc


    def hydrate_blinded_batch_id(self, bundle):
        if bundle.data.get("blinded_batch_id", "") == u"EMPTY_ID":
            uox_id = generate_uox_id()
            bundle.data["blinded_batch_id"] = uox_id
            bundle.obj.blinded_batch_id = uox_id
        return bundle


    def dehydrate_properties(self, bundle):
        archived = bundle.obj.properties.get("archived", "false")
        value = bundle.obj.properties
        value["archived"] = archived
        return value


    def get_project_specific_data(self, request, queries, pids, sorts, textsearch):
        to_return = []
        for pid in pids:
            #We have removed the schema from the index for data efficiency so
            #We need to pull it out of the project resource
            pr = ChemregProjectResource()
            obj = pr._meta.queryset.get(pk=pid)
            bundle = pr.build_bundle(obj=obj, request=request)
            bundle = pr.full_dehydrate(bundle)
            
            proj_data = self._meta.serializer.to_simple(bundle.data, {})
            project_data_fields_json = proj_data["custom_field_config"]["project_data_fields"]

            table_schema_for_project = reformat_project_data_fields_as_table_schema( "export", project_data_fields_json)

            proj_index = elasticsearch_client.get_project_index_name(pid)
            es_data = elasticsearch_client.get_list_data_elasticsearch(queries,
                    proj_index,
                    sorts=sorts, 
                    textsearch=textsearch,
                    offset=0,
                    limit=10000 )
            standardised = self.prepare_es_hits(es_data)

            standardised["name"] = proj_data["name"]
            standardised["id"] = pid
            standardised["schema"] = table_schema_for_project
            to_return.append(standardised)
        query_reps = [{"fieldn": "Project", "pick_from_list": ",".join([str(ret["name"]) for ret in to_return]) }]
        for i, q in enumerate(queries):
            qrep = {qtype["value"] : u"" for qtype in settings.CBH_QUERY_TYPES}
            qrep["fieldn"] = ""
            for key, value in q.items():
                if key == "pick_from_list":
                    qrep[key] = ", ".join(value)
                elif key == "field_path":
                    if value:
                        knownBy = settings.TABULAR_DATA_SETTINGS["schema"].get(value,{}).get("knownBy", "")
                        if not knownBy:
                            knownBy = str(value.split(".")[-1])
                        qrep["fieldn"] = "".join([l for l in knownBy])
                else:
                    qrep[key] = str(value)

            qrep["id"] = i+1
            query_reps.append(qrep)

        schem = [{"data": "id", "knownBy" : "Row"}, {"data": "fieldn", "knownBy" : "Col"}] + [{"data": qtype["value"], "knownBy" : qtype["name"]} for qtype in settings.CBH_QUERY_TYPES]
        query_summary = {
                "name" : "_Query Used",
                "objects" : query_reps,
                "schema" : schem
                }
        to_return.append(query_summary)
        return to_return



    def get_detail(self, request, **kwargs):
        """
        Returns a single serialized resource.
        Calls ``cached_obj_get/obj_get`` to provide the data, then handles that result
        set and serializes it.
        Should return a HttpResponse (200 OK).
        """
        basic_bundle = self.build_bundle(request=request)

        kwargs = self.remove_api_resource_names(kwargs)
        allowed_pids = self._meta.authorization.project_ids(request)
        concatenated_indices = elasticsearch_client.get_list_of_indicies(allowed_pids)
        data = elasticsearch_client.get_detail_data_elasticsearch(concatenated_indices, kwargs["pk"])   
        

        requested_pid = data["projectfull"]["id"]

        if requested_pid not in allowed_pids:
            raise Unauthorized("No permissions for requested project") 

        bundle = self.alter_detail_data_to_serialize(request, data)
        return self.create_response(request, bundle)


    def get_list(self, request, **kwargs):
        """
        Returns a serialized list of resources.
        Calls ``obj_get_list`` to provide the data, then handles that result
        set and serializes it.
        Should return a HttpResponse (200 OK).
        """
        # TODO: Uncached for now. Invalidation that works for everyone may be
        #       impossible.

        base_bundle = self.build_bundle(request=request)
        pids = request.GET.get("pids", "")
        project_ids = []
        if pids:

            project_ids = [int(pid) for pid in pids.split(",")]
        
        allowed_pids = set(self._meta.authorization.project_ids(request))

        
        for requested_pid in project_ids:
            if requested_pid not in allowed_pids:
                raise Unauthorized("No permissions for requested project") 



        if len(project_ids) == 0:
            project_ids = allowed_pids

        queries = json.loads(b64decode(request.GET.get("encoded_query", EMPTY_ARRAY_B64)))

        extra_queries = kwargs.get("extra_queries", False)
        if extra_queries:
            #extra queries can be added via kwargs
            queries += extra_queries

        #Search for whether this item is archived or not ("archived is indexed as a string")
        archived = request.GET.get("archived", "false")
        queries.append({"query_type": "phrase", "field_path": "properties.archived", "phrase": archived})

        sorts = json.loads(b64decode(request.GET.get("encoded_sorts", EMPTY_ARRAY_B64)))
        if len(sorts) == 0:
            sorts = [{"field_path":"id","sort_direction":"desc"}]
        textsearch = b64decode(request.GET.get("textsearch", ""))
        limit = request.GET.get("limit", 10)
        offset = request.GET.get("offset", 0)
        autocomplete = request.GET.get("autocomplete", "")
        autocomplete_field_path = request.GET.get("autocomplete_field_path", "")
        autocomplete_size = request.GET.get("autocomplete_size", settings.MAX_AUTOCOMPLETE_SIZE)

        concatenated_indices = elasticsearch_client.get_list_of_indicies(project_ids)

        if request.GET.get("format", None) != "sdf" and request.GET.get("format", None) != "xlsx":
            data = elasticsearch_client.get_list_data_elasticsearch(queries,
                concatenated_indices,
                sorts=sorts, 
                offset=offset, 
                limit=limit, 
                textsearch=textsearch, 
                autocomplete=autocomplete,
                autocomplete_field_path=autocomplete_field_path,
                autocomplete_size=autocomplete_size )
        else:
            #Limit , offset and autocomplete have no effect for a project export
            data = self.get_project_specific_data(request,  queries, project_ids, sorts, textsearch)
            
            return self.create_response(request, data) 

        if autocomplete_field_path:
            
            bucks = data["aggregations"]["filtered_field_path"]["field_path_terms"]["buckets"]
            bundledata = {"items" : bucks,
                            "autocomplete" : autocomplete,
                            "unique_count" : data["aggregations"]["filtered_field_path"]["unique_count"]["value"]}
        else:
            #This is just a standard request for data
            bundledata = self.prepare_es_hits(data)

            bundledata = self.alter_list_data_to_serialize(request, bundledata)
        return self.create_response(request, bundledata) 

    def prepare_es_hits(self, hits):
        return {"objects": 
                            [hit["_source"] for hit in hits["hits"]["hits"]],
                            "meta" : {"total_count" : hits["hits"]["total"]}
                            }

    class Meta:
        authorization = ProjectAuthorization()
        queryset = CBHCompoundBatch.objects.all()
        
        include_resource_uri = True
        serializer = CBHCompoundBatchSerializer()
        always_return_data = True







class CBHCompoundBatchResource(BaseCBHCompoundBatchResource):

    class Meta(BaseCBHCompoundBatchResource.Meta):
        resource_name = 'cbh_compound_batches_v2'


    
    def get_list(self, request, **kwargs):
        """Request data where the project is not of the saved search type"""
        return super(CBHCompoundBatchResource, self).get_list(request, extra_queries=[{"query_type": "pick_from_list", "field_path" : "projectfull.project_type.saved_search_project_type", "pick_from_list" : ["false"] }])



class CBHSavedSearchResource(BaseCBHCompoundBatchResource):
    class Meta(BaseCBHCompoundBatchResource.Meta):
        resource_name = 'cbh_saved_search'


    def get_list(self, request, **kwargs):
        """Request data where the project is of the saved search type"""

        return super(CBHSavedSearchResource, self).get_list(request, extra_queries=[{"query_type": "pick_from_list", "field_path" : "projectfull.project_type.saved_search_project_type", "pick_from_list" : ["true"] }])




class IndexingCBHCompoundBatchResource(BaseCBHCompoundBatchResource):

    def prepare_fields_for_index(self, batch_dicts):
        """Fields are stored as strings in hstore so convert the list or object fields back from JSON"""
        for bd in batch_dicts:
            for field in bd["projectfull"]["custom_field_config"]["project_data_fields"]:
                
                if field["edit_schema"]["properties"][field["name"]]["type"] == "object":
                    if bd["custom_fields"].get(field["name"], False):
                        bd["custom_fields"][field["name"]] = json.loads(bd["custom_fields"][field["name"]])
                    else:
                        bd["custom_fields"][field["name"]] = field["edit_schema"]["properties"][field["name"]]["default"]
                if field["edit_schema"]["properties"][field["name"]]["type"] == "array":
                    if bd["custom_fields"].get(field["name"], False):
                        bd["custom_fields"][field["name"]] = json.loads(bd["custom_fields"][field["name"]])
                    else:
                        bd["custom_fields"][field["name"]] = []
        return batch_dicts


    def index_batch_list(self, request, batch_list):
        #retrieve some objects as json
        bundles = [
            self.full_dehydrate(self.build_bundle(obj=obj, request=request), for_list=True)
            for obj in batch_list
        ]

        #retrieve schemas which tell the elasticsearch request which fields to index for each object (we avoid deserializing a single custom field config more than once)
        #Now make the schema list parallel to the batches list
        batch_dicts = [self.Meta.serializer.to_simple(bun, {}) for bun in bundles]

        batch_dicts = self.prepare_fields_for_index(batch_dicts)

        project_data_field_sets = [batch_dict["projectfull"]["custom_field_config"].pop("project_data_fields") for batch_dict in batch_dicts]
        schemas = [reformat_project_data_fields_as_table_schema( "indexing", pdfs) for pdfs in project_data_field_sets]
        index_names = []

        for batch in batch_dicts:
            batch["projectfull"]["custom_field_config"] = batch["projectfull"]["custom_field_config"]["resource_uri"]
            index_name = elasticsearch_client.get_project_index_name(batch["projectfull"]["id"])
            index_names.append(index_name)
        
        batch_dicts = elasticsearch_client.index_dataset(batch_dicts, schemas, index_names)
            

    def reindex_elasticsearch(self, request, **kwargs):
        desired_format = self.determine_format(request)
        batches = self.get_object_list(request)
        # we only want to store certain fields in the search index
        from django.core.paginator import Paginator
        paginator = Paginator(batches, 1000) # chunks of 1000

        for page in range(1, paginator.num_pages +1):
            bs = paginator.page(page).object_list
            self.index_batch_list(request, bs)
            
            # here you can do what you want with the row
            
            print "done page %d of %d" % (page ,paginator.num_pages)

       
        return HttpResponse(content="test", )
        
    class Meta(BaseCBHCompoundBatchResource.Meta):
        pass


def index_batches_in_new_index(batches):
    """Temporary function to index all the data"""
    request = HttpRequest()
    IndexingCBHCompoundBatchResource().index_batch_list(request, batches)




##Multiple batch creation


class CBHMultipleBatchUploadResource(ModelResource):
    pass

