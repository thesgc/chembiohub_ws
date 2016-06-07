"""
nearly deprecated old inteface to the compounds API
"""
from tastypie.resources import ALL
from tastypie.resources import ALL_WITH_RELATIONS
from tastypie.resources import ModelResource
from django.conf import settings
from django.conf.urls import *
from django.core.exceptions import ObjectDoesNotExist
from tastypie.authorization import Authorization
from tastypie import http
from django.http import HttpResponse
from collections import OrderedDict
from pybel import readfile, readstring
import re
import shortuuid
import copy
from cbh_utils import elasticsearch_client
import math
from tastypie.exceptions import ImmediateHttpResponse

try:
    from rdkit import Chem
except ImportError:
    Chem = None
from tastypie.exceptions import BadRequest
from django_q.tasks import async_iter, result, async
try:
    from chembl_compatibility.models import MoleculeDictionary
    from chembl_compatibility.models import CompoundMols
except ImportError:
    from chembl_core_model.models import MoleculeDictionary
    from chembl_core_model.models import CompoundMols
try:
    WS_DEBUG = settings.WS_DEBUG
except AttributeError:
    WS_DEBUG = False
from cbh_core_api.authorization import ProjectAuthorization
from cbh_core_api.resources import ChemregProjectResource, ChemRegCustomFieldConfigResource, NoCustomFieldsChemregProjectResource
from tastypie.utils import dict_strip_unicode_keys
from tastypie.serializers import Serializer
from tastypie import fields
from cbh_chembl_model_extension.models import CBHCompoundBatch, CBHCompoundMultipleBatch, generate_structure_and_dictionary
from cbh_core_model.models import Project, PinnedCustomField
from tastypie.authentication import SessionAuthentication
import json
from tastypie.paginator import Paginator
from cbh_core_model.models import CBHFlowFile
import pandas as pd
from django.db.models import Max, Q
from tastypie.authentication import SessionAuthentication
import chemdraw_reaction
from django.contrib.auth import get_user_model
from django.http import HttpRequest
from rdkit.Chem.AllChem import Compute2DCoords
from django.db.models import Prefetch
import dateutil.parser
from cbh_utils.parser import parse_pandas_record, parse_sdf_record, apply_json_patch, get_uncurated_fields_from_file
# from tastypie.utils.mime import build_content_type
from cbh_core_api.resources import SimpleResourceURIField, UserResource, UserHydrate
import time
from cbh_chem_api.tasks import  get_batch_from_sdf_chunks, get_batch_from_xls_row, process_file_request, validate_multi_batch
import itertools
from tastypie.http import  HttpConflict


def build_content_type(format, encoding='utf-8'):
    """
    Appends character encoding to the provided format if not already present.
    """
    if 'charset' in format:
        return format
    return "%s; charset=%s" % (format, encoding)





class CBHCompoundUploadResource(ModelResource):
    """Old CBHCompoundBatch resource now only used for bulk data uploads"""
    project = SimpleResourceURIField(
        NoCustomFieldsChemregProjectResource, 'project_id', blank=False, null=False, help_text="")
    creator = SimpleResourceURIField(UserResource, 'created_by_id', null=True, readonly=True)

    class Meta:
        
        always_return_data = True
        prefix = "related_molregno"
        queryset = CBHCompoundBatch.objects.all()
        resource_name = 'cbh_compound_batches'
        authorization = ProjectAuthorization()
        include_resource_uri = False
        always_return_data = True
        serializer = Serializer()
        allowed_methods = ['get', 'post', 'put', 'patch']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        paginator_class = Paginator
        # elasticsearch config items



    def alter_deserialized_list_data(self, request, deserialized):
        """Weird way we added project to input data, not reccomended, now deprecated apart from for upload"""
        proj = Project.objects.get(project_key=deserialized["project_key"])
        deserialized["project"] = proj
        return deserialized

    def alter_deserialized_detail_data(self, request, deserialized):
        '''A project may be necessary for create statements'''
        if deserialized["project_key"]:
            proj = Project.objects.get(project_key=deserialized["project_key"])
            deserialized["project"] = proj
        return deserialized

    def full_hydrate(self, bundle):
        '''As the object is created we run the validate code on it'''
        bundle = super(CBHCompoundUploadResource, self).full_hydrate(bundle)
        if not bundle.obj.id:
            bundle.obj.created_by = bundle.request.user.username
            bundle.obj.created_by_id = bundle.request.user.id
            # bundle.obj.validate()
            # self.match_list_to_moleculedictionaries(bundle.obj,bundle.data["project"] )
        return bundle

    def prepend_urls(self):
        """Add lots of extra methods to the API"""
        return [
            url(r"^(?P<resource_name>%s)/delete_index/?$" % self._meta.resource_name,
                self.wrap_view('delete_index'), name="delete_index"),
            url(r"^(?P<resource_name>%s)/update_temp_batches/?$" % self._meta.resource_name,
                self.wrap_view('update_temp_batches'), name="update_temp_batches"),
            url(r"^(?P<resource_name>%s)/get_part_processed_multiple_batch/?$" % self._meta.resource_name,
                self.wrap_view('get_part_processed_multiple_batch'), name="api_get_part_processed_multiple_batch"),

            url(r"^(?P<resource_name>%s)/validate_list/?$" % self._meta.resource_name,
                self.wrap_view('post_validate_list'), name="api_validate_compound_list"),
            url(r"^(?P<resource_name>%s)/multi_batch_save/?$" % self._meta.resource_name,
                self.wrap_view('multi_batch_save'), name="multi_batch_save"),
            url(r"^(?P<resource_name>%s)/multi_batch_custom_fields/?$" % self._meta.resource_name,
                self.wrap_view('multi_batch_custom_fields'), name="multi_batch_custom_fields"),
            url(r"^(?P<resource_name>%s)/validate_files/?$" % self._meta.resource_name,
                self.wrap_view('post_validate_files'), name="api_compound_validate_files"),]


  
    def patch_dict(self, dictdata, headers):
        """Apply a JSON patch to the data in a particular column"""
        for header in headers:
            json_patches = copy.copy(header.get("operations", False))
            if json_patches:
                apply_json_patch(dictdata, json_patches)

    def multi_batch_save(self, request, **kwargs):
        """Save the data which has been cached in Elasticsearch"""
        deserialized = self.deserialize(request, request.body, format=request.META.get(
            'CONTENT_TYPE', 'application/json'))
        session_key = request.COOKIES[settings.SESSION_COOKIE_NAME]
        deserialized = self.alter_deserialized_detail_data(
            request, deserialized)

        bundle = self.build_bundle(
            data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(
                self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(
                self.get_object_list(bundle.request), bundle)
        
        id = bundle.data["multiplebatch"]
        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        creator_user = request.user
        if not bundle.data.get("task_id_for_save", None):

            mb.created_by = creator_user.username

            bundle.data["task_id_for_save"] = async('cbh_chem_api.tasks.save_multiple_batch',  mb, creator_user, session_key)
           
            
        
        res = result(bundle.data["task_id_for_save"], wait=300)
        if res is True:
            return self.create_response(request, bundle, response_class=http.HttpCreated)
        return self.create_response(request, bundle, response_class=http.HttpAccepted)
        


    def after_save_and_index_hook(self, request, multi_batch_id, project_id):
        """Hook used to perform operations on data that has been saved"""
        pass

    def alter_batch_data_after_save(self, batch_list, python_file_object, multi_batch):
        """Actually edit the data just after it has been saved"""
        pass


    def delete_index(self, request, **kwargs):
        """Delete the index that was created for a multiple batch"""
        deserialized = self.deserialize(request, request.body, format=request.META.get(
            'CONTENT_TYPE', 'application/json'))

        deserialized = self.alter_deserialized_detail_data(
            request, deserialized)
        session_key = request.COOKIES[settings.SESSION_COOKIE_NAME]
        bundle = self.build_bundle(
            data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(
                self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(
                self.get_object_list(bundle.request), bundle)
        id = bundle.data["multiplebatch"]
        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        elasticsearch_client.delete_index(
            elasticsearch_client.get_temp_index_name(session_key, mb.id))
        return self.create_response(request, bundle, response_class=http.HttpAccepted)

   

    def update_temp_batches(self, request, **kwargs):
        '''Update a set of molecules into elasticsearch (used in ChemBio Hub to set the action field to ignore or new batch)'''
        deserialized = self.deserialize(request, request.body, format=request.META.get(
            'CONTENT_TYPE', 'application/json'))

        deserialized = self.alter_deserialized_detail_data(
            request, deserialized)
        bundle = self.build_bundle(
            data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(
                self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(
                self.get_object_list(bundle.request), bundle)
        multi_batch_id = bundle.data["multiplebatch"]
        es_ready_updates = bundle.data["objects"]

        index_name = elasticsearch_client.get_temp_index_name(
            request.COOKIES[settings.SESSION_COOKIE_NAME], multi_batch_id)

        elasticsearch_client.create_temporary_index(
            es_ready_updates, index_name)
        elasticsearch_client.get_action_totals(index_name, bundle.data)
        return self.create_response(request, bundle, response_class=http.HttpAccepted)

    def multi_batch_custom_fields(self, request, **kwargs):
        '''change the structure column for an excel file'''
        deserialized = self.deserialize(request, request.body, format=request.META.get(
            'CONTENT_TYPE', 'application/json'))

        deserialized = self.alter_deserialized_detail_data(
            request, deserialized)
        bundle = self.build_bundle(
            data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(
                self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(
                self.get_object_list(bundle.request), bundle)
        id = bundle.data["multiplebatch"]
        headers = bundle.data["headers"]
        # structure_col = bundle.data.get("structure_col", None)
        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        processSmiles = False
        # if structure_col and structure_col != mb.uploaded_data.get("structure_col", ""):
        #     processSmiles =  True
        index_name = elasticsearch_client.get_temp_index_name(request, mb.id)
        elasticsearch_client.get_action_totals(index_name, bundle.data)
        mb.uploaded_data = bundle.data
        mb.save()
        return self.create_response(request, bundle, response_class=http.HttpAccepted)

    

        

    def get_part_processed_multiple_batch(self, request, **kwargs):
        """
        Get the part processed data from elasticsearch and the stats about the
        multiple batch
        """
        # TODO: Uncached for now. Invalidation that works for everyone may be
        #       impossible.
        bundle = self.build_bundle(request=request)
        session_key = request.COOKIES[settings.SESSION_COOKIE_NAME]
        # self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        if(kwargs.get("multi_batch", None)):
            mb = kwargs.get("multi_batch")
            id = mb.id
        else:
            id = request.GET.get("current_batch")
            mb = CBHCompoundMultipleBatch.objects.get(pk=id)

        if not mb.uploaded_data:
            #The uploaded data field will be set once the data is fully processed
            raise ImmediateHttpResponse(HttpConflict("data_not_yet_ready"))


        to_be_serialized = mb.uploaded_data
        to_be_serialized = self.get_cached_temporary_batch_data(
            id, request.GET, session_key, bundledata=to_be_serialized)
        index_name = elasticsearch_client.get_temp_index_name(session_key, id)
        elasticsearch_client.get_action_totals(index_name, to_be_serialized)
        return self.create_response(request, to_be_serialized)

    def post_validate_list(self, request, **kwargs):
        """
        When a list of SMILES patterns are submitted, save them to elasticsearhc and
        call the normal validate multi batch function to give the normal statistics page
        """
        session_key = request.COOKIES[settings.SESSION_COOKIE_NAME]
        deserialized = self.deserialize(request, request.body, format=request.META.get(
            'CONTENT_TYPE', 'application/json'))

        deserialized = self.alter_deserialized_detail_data(
            request, deserialized)
        bundle = self.build_bundle(
            data=dict_strip_unicode_keys(deserialized), request=request)
        if bundle.obj.pk:
            self.authorized_update_detail(
                self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(
                self.get_object_list(bundle.request), bundle)
        smilesdata = bundle.data.get("smilesdata", "")
        objects = [smi.strip() for smi in smilesdata.splitlines(True)]
        batches = []
        # first assume smiles
        allmols = [(obj, Chem.MolFromSmiles(str(obj))) for obj in objects]
        # Next test for inchi
        multiple_batch = CBHCompoundMultipleBatch.objects.create(
            project=bundle.data["project"], batch_count=len(allmols))
        for index, m in enumerate(allmols):
            if m[1] is None:
                inchimol = Chem.MolFromInchi(
                    str(m[0].encode('ascii', 'ignore')))
                if inchimol is not None:
                    allmols[index] = (Chem.MolToSmiles(inchimol), inchimol)
        for m in allmols:
            if(m[1]):
                Compute2DCoords(m[1])
        batches = []
        
        for mol2 in allmols:
            if mol2[1]:
                b = CBHCompoundBatch.objects.from_rd_mol(
                    mol2[1], smiles=mol2[0], project=bundle.data["project"], reDraw=True)
            else:
                b = CBHCompoundBatch.objects.blinded(
                    project=bundle.data["project"])
                b.warnings["smilesParseError"] = "true"
                b.properties["action"] = "Ignore"
                b.original_smiles = mol2[0]
            b.multiple_batch_id = multiple_batch.pk
            b.created_by = bundle.request.user.username
            b.created_by_id = bundle.request.user.id
            batches.append(b)

        bundle.data["current_batch"] = multiple_batch.pk
        bundle.data["headers"] = []
        validate_multi_batch(self, multiple_batch, bundle.data, session_key, batches)
        return self.create_response(request, bundle, response_class=http.HttpAccepted)

    
    def preprocess_sdf_file(self, file_obj, request):
        """Read the input SDF file and add some extra information to it"""
        pass



    def post_validate_files(self, request, **kwargs):
        """Receive a CBHFlowFile ID which points at an uploaded SDF, XLSX or CDX file
        Perform validation on the file's contents and then send the resultant elasticsearch index
        to the validate mult batch method
        More about data import information can be found in the wiki
        """
        automapped_structure = False
        deserialized = self.deserialize(request, request.body, format=request.META.get(
            'CONTENT_TYPE', 'application/json'))

        deserialized = self.alter_deserialized_detail_data(
            request, deserialized)
        bundle = self.build_bundle(
            data=dict_strip_unicode_keys(deserialized), request=request)
        self.authorized_create_detail(
            self.get_object_list(bundle.request), bundle)
        file_name = bundle.data['file_name']
        session_key = request.COOKIES[settings.SESSION_COOKIE_NAME]
        correct_file = CBHFlowFile.objects.get(
            identifier="%s-%s" % (session_key, file_name))

        if(correct_file.extension in (".xls", ".xlsx")):
            # we need to know which column contains structural info - this needs to be defined on the mapping page and passed here
            # read in the specified structural column

            df = None
            try:
                df = pd.read_excel(correct_file.file)
            except IndexError:
                raise BadRequest("no_headers")
        if(correct_file.extension not in (".xls", ".xlsx", ".sdf", ".cdxml")):
            raise BadRequest("file_format_error")
        print correct_file.file.name
        bundledata = bundle.data
        creator_user = request.user
        cfr = ChemRegCustomFieldConfigResource()
        jsondata = json.loads(cfr.get_detail(request,
            pk=bundledata["project"].custom_field_config_id).content)
        schemaform = [field["edit_form"]["form"][0] for field in jsondata["project_data_fields"]]
        if (correct_file.extension == ".sdf"):
            # read in the filepalter_batch_data_after_save(self, batch_list, python_file_obj, request, multi_batch):
            self.preprocess_sdf_file(correct_file.file, request)
        multiple_batch = CBHCompoundMultipleBatch.objects.create(
            project=bundledata["project"],
            uploaded_file=correct_file
        )
        bundle.data["current_batch"] = multiple_batch.pk
        bundle.data["multiplebatch"] = multiple_batch.pk
        id = async("cbh_chem_api.tasks.process_file_request", 
                                                               multiple_batch,
                                                               bundledata, 
                                                               creator_user,
                                                               schemaform,
                                                               correct_file,
                                                               session_key)

        # res = result(id, wait=20000)
        
        bundle.data = bundledata
        return self.create_response(request, bundle, response_class=http.HttpAccepted)


    def dehydrate(self, bundle):
        """Tidy up the data adding timestamp etc, now mostly deprecated"""
        # try:

        mynames = [ "uncurated_fields", "warnings", "properties", "custom_fields",]
        for name in mynames:
            bundle.data[name] = json.loads(bundle.data[name])
            
       
        return bundle

    def batches_to_es_ready(self, batches, non_chem_data_only=None):
        """Convert the data to be ready to submit to elasticsearch"""
        batch_dicts = []
        index = 1
        for batch in batches:
            if batch:
                if(not batch.id):
                    batch.id = index
                bun = self.build_bundle(obj=batch, request=HttpRequest())
                bun = self.full_dehydrate(bun, for_list=True)
                if bun.data["project"] == '':
                    bun.data[
                        "project"] = '/%s/cbh_projects/%d' % (settings.WEBSERVICES_NAME, bun.obj.project_id)
                if non_chem_data_only:
                    ready = bun.data
                else:
                    ready = bun.data

                batch_dicts.append(ready)
            else:
                # preserve the line number of the batch that could not be
                # processed
                batch_dicts.append(
                    {"id": index,
                     "warnings": {"parseerror": "true"
                                  },
                     "properties": {"action": "Ignore"},
                     "project": str(self.project)
                     }
                )
            index += 1
        return batch_dicts

    def set_cached_temporary_batches(self, batches, multi_batch_id, session_key):
        """Index the new data when a new bulk upload is done"""
        batch_dicts = self.batches_to_es_ready(batches)
        index_name = elasticsearch_client.get_temp_index_name(
            session_key, multi_batch_id)
        elasticsearch_client.create_temporary_index(
            batch_dicts,  index_name)
        # Now get rid of my ES preparation again

    def get_cached_temporary_batches(self, bundles, request, bundledata={}):
        """Request the imported data which has been cached in Elasticsearch"""
        fakeobj = CBHCompoundBatch(project=bundledata["project"], id=10)
        bun = self.build_bundle(obj=fakeobj)
        bun = self.full_dehydrate(bun)
        proj = bun.data["project"]
        data = []
        for datum in bundles["objects"]:
        
            datum = self.build_bundle(data=datum, request=request)

            datum.data["project"] = proj
            datum = self.full_hydrate(datum)
            #Quick hack to fix saving of sutom fields when validating a sketch
            if bundledata.get("type") == "sketch" and len(bundles["objects"]) == 1:
                datum.data["custom_fields"] = bundledata["custom_fields"]
                datum.obj.custom_fields = bundledata["custom_fields"]
            data.append(datum)
        bundles["objects"] = data
        return bundles

   
    def get_cached_temporary_batch_data(self, multi_batch_id, get_data, session_key, bundledata={}):
        """make the batch data into models so it can be serialized properly"""
        es_request = {
            "from": get_data.get("offset", 0),
            "size": get_data.get("limit", 50),
            "filter": json.loads(get_data.get("query", '{ "match_all" : {}}')),
            "sort": json.loads(get_data.get("sorts", '[{"id": {"order": "asc"}}]'))
        }
        index = elasticsearch_client.get_temp_index_name(
            session_key, multi_batch_id)
        bundledata = elasticsearch_client.get_from_temp_index(index, es_request, bundledata)
        
        return bundledata



def deepgetattr(obj, attr, ex):
    """Recurses through an attribute chain to get the ultimate value."""
    try:
        return reduce(getattr, attr.split('.'), obj)
    except:
        return ex


    



