from chembl_business_model.models.compounds import CompoundProperties, MoleculeDictionary
from cbh_chembl_model_extension.models import CBHCompoundBatch
from tastypie.resources import ALL
from tastypie.resources import ALL_WITH_RELATIONS
from tastypie.resources import ModelResource
from tastypie.serializers import Serializer
from django.conf import settings
from django.conf.urls import *
from django.core.exceptions import ObjectDoesNotExist
from tastypie.authorization import Authorization
from tastypie import fields
from cbh_chem_api.projects import ChemregProjectResource, ChemRegCustomFieldConfigResource, NoCustomFieldsChemregProjectResource
from cbh_core_api.resources import SimpleResourceURIField, UserResource, UserHydrate, CBHNoSerializedDictField
from cbh_utils import elasticsearch_client
import json 
from django.http import HttpResponse





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
         NoCustomFieldsChemregProjectResource, 'project', blank=False, null=False, full=True, readonly=True)
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









class IndexingCBHCompoundBatchResource(BaseCBHCompoundBatchResource):

    def fix_data_types_for_index(self, value):
        """Elasticsearch will not index dictionaries"""
        if not value:
            return "__EMPTY"
        if isinstance(value, basestring):
            return value
        if type(value) is dict:
            return json.dumps(value)
        if type(value) is list:
            return [self.fix_data_types_for_index(v) for v in value]
        return unicode(value)

    def yield_indexed_fields(self, bundle):
        for field in bundle.data["custom_field_config_full"]["project_data_fields"]:
            field_value = bundle.data["custom_fields"].get(field["name"], None)
            yield {"name": field["name"], "value": self.fix_data_types_for_index(field_value)}

    def alter_list_data_to_serialize(self, request, data):
        """Here we get a cached copy of the custom field config and use it to ensure each field is properly indexed"""
        custom_fields_to_retrieve = set()
        for bun in data[self._meta.collection_name]:
            custom_fields_to_retrieve.add(bun.obj.project.custom_field_config_id)
        cfcr = ChemRegCustomFieldConfigResource()

        bundles = [
            cfcr.full_dehydrate(cfcr.build_bundle(obj=obj, request=request), for_list=True)
            for obj in cfcr.Meta.queryset.filter(id__in=custom_fields_to_retrieve)
        ]

        simple_dicts = [cfcr.Meta.serializer.to_simple(bun, {}) for bun in bundles]


        cfc_lookup = { cfc["resource_uri"] : cfc for cfc in simple_dicts }
        for bun in data[self._meta.collection_name]:
            elasticsearch_client.prepare_dataset_for_indexing
            bun.data["custom_field_config_full"] = cfc_lookup[bun.data["projectfull"].data["custom_field_config"]]
            bun.data["indexed_fields"] = list(self.yield_indexed_fields(bun))
            del bun.data["custom_field_config_full"]



    def reindex_elasticsearch(self, request, **kwargs):
        desired_format = self.determine_format(request)
        batches = self.get_object_list(request)
        # we only want to store certain fields in the search index
        from django.core.paginator import Paginator
        paginator = Paginator(batches, 1000) # chunks of 1000

        for page in range(1, paginator.num_pages +1):
            bs = paginator.page(page).object_list
            bundles = [
                self.full_dehydrate(self.build_bundle(obj=obj, request=request), for_list=True)
                for obj in bs
            ]
            self.alter_list_data_to_serialize(request, {"objects" : bundles} )
            batch_dicts = [self.Meta.serializer.to_simple(bun, {}) for bun in bundles]

            # reindex compound data
            from pprint import pprint
            pprint(batch_dicts)
            index_name = elasticsearch_client.get_main_index_name()
            es_reindex = elasticsearch_client.create_temporary_index(
                batch_dicts, request, index_name)
            if es_reindex.get("errors"):
                print "ERRORS"
                print json.dumps(es_reindex)
            # here you can do what you want with the row
            
            print "done page %d of %d" % (page ,paginator.num_pages)

       
        return HttpResponse(content="test", )
        
    class Meta(BaseCBHCompoundBatchResource.Meta):
        pass






