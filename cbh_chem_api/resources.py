from chembl_business_model.models.compounds import CompoundProperties, MoleculeDictionary
from cbh_chembl_model_extension.models import CBHCompoundBatch
from tastypie.resources import ALL
from tastypie.resources import ALL_WITH_RELATIONS
from tastypie.resources import ModelResource, Resource
from tastypie.serializers import Serializer
from django.conf import settings
from django.conf.urls import *
from django.core.exceptions import ObjectDoesNotExist
from tastypie.authorization import Authorization
from tastypie import fields
from cbh_core_api.resources import SimpleResourceURIField, UserResource, UserHydrate, CBHNoSerializedDictField,  ChemregProjectResource, ChemRegCustomFieldConfigResource, NoCustomFieldsChemregProjectResource
from cbh_utils import elasticsearch_client
import json 
from django.http import HttpResponse, HttpRequest


class CompoundPropertiesResource(ModelResource):
    class Meta:
        queryset = CompoundProperties.objects.all()
        fields = ["alogp", "full_molformula", "full_mwt"]

class MoleculeDictionaryResource(ModelResource):
    compoundproperties = fields.ForeignKey(CompoundPropertiesResource, 'compoundproperties',  null=True, readonly=True, full=True)

    class Meta:
        queryset = MoleculeDictionary.objects.all()
        fields = ["compoundproperties"]


class BaseCBHCompoundBatchResource(ModelResource):
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

    class Meta:
        queryset = CBHCompoundBatch.objects.all()
        resource_name = 'cbh_compound_batches'
        include_resource_uri = True
        serializer = Serializer()


class CBHCompoundBatchSearchResource(Resource):

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

        queries = json.loads(request.GET.get("encoded_query", "[]"))
        limit = request.GET.get("limit", 10)
        offset = request.GET.get("offset", 0)
        index = elasticsearch_client.get_main_index_name()
        data = elasticsearch_client.get_list_data_elasticsearch(queries,index, offset=offset, limit=limit )
        bundledata = {"objects": 
                        [hit["_source"] for hit in data["hits"]["hits"]],
                        "meta" : {"totalCount" : data["hits"]["total"]}
                        }
        bundledata = self.alter_list_data_to_serialize(request, bundledata)
        return self.create_response(request, bundledata) 


    class Meta:
        resource_name = "cbh_compound_batches_search"



class IndexingCBHCompoundBatchResource(BaseCBHCompoundBatchResource):

    def reformat_project_data_fields_as_table_schema(self, table_schema_type, project_data_fields_json):
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

    def index_batch_list(self, request, batch_list):

            #retrieve some objects as json
        bundles = [
            self.full_dehydrate(self.build_bundle(obj=obj, request=request), for_list=True)
            for obj in batch_list
        ]

        #retrieve schemas which tell the elasticsearch request which fields to index for each object (we avoid deserializing a single custom field config more than once)
        #Now make the schema list parallel to the batches list
        batch_dicts = [self.Meta.serializer.to_simple(bun, {}) for bun in bundles]

        project_data_field_sets = [batch_dict["projectfull"]["custom_field_config"].pop("project_data_fields") for batch_dict in batch_dicts]

        schemas = [self.reformat_project_data_fields_as_table_schema( "indexing", pdfs) for pdfs in project_data_field_sets]
        for batch in batch_dicts:
            batch["projectfull"]["custom_field_config"] = batch["projectfull"]["custom_field_config"]["resource_uri"]
        index_name = elasticsearch_client.get_main_index_name()
        batch_dicts = elasticsearch_client.index_dataset(index_name, batch_dicts, schemas)
            

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


