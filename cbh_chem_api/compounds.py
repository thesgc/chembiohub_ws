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
import elasticsearch_client
from difflib import SequenceMatcher as fuzzymatch
import re
try:
    from rdkit import Chem
except ImportError:
    Chem = None
from tastypie.exceptions import BadRequest

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
from cbh_chem_api.serializers import CBHCompoundBatchSerializer, CBHCompoundBatchElasticSearchSerializer, get_key_from_field_name
from tastypie.utils import dict_strip_unicode_keys
from tastypie.serializers import Serializer
from tastypie import fields
from cbh_chembl_model_extension.models import CBHCompoundBatch, CBHCompoundMultipleBatch
from cbh_core_model.models import Project, PinnedCustomField
from tastypie.authentication import SessionAuthentication
import json
from tastypie.paginator import Paginator
from flowjs.models import FlowFile
import pandas as pd
from django.db.models import Max, Q
from tastypie.serializers import Serializer
from tastypie.authentication import SessionAuthentication
import chemdraw_reaction
from django.contrib.auth import get_user_model
from rdkit.Chem.AllChem import Compute2DCoords
from django.db.models import Prefetch
import dateutil.parser
from cbh_utils.parser import parse_pandas_record, parse_sdf_record, apply_json_patch, get_uncurated_fields_from_file
# from tastypie.utils.mime import build_content_type
from cbh_core_api.resources import SimpleResourceURIField, UserResource, UserHydrate
import time
from cbh_chem_api.resources import index_batches_in_new_index

def build_content_type(format, encoding='utf-8'):
    """
    Appends character encoding to the provided format if not already present.
    """
    if 'charset' in format:
        return format
    return "%s; charset=%s" % (format, encoding)




class CBHCompoundBatchResource(ModelResource):
    project = fields.ForeignKey(
        ChemregProjectResource, 'project', blank=False, null=False)
    creator = SimpleResourceURIField(UserResource, 'created_by_id', null=True, readonly=True)
    projectfull = fields.ForeignKey(
         NoCustomFieldsChemregProjectResource, 'project', blank=False, null=False, full=True, readonly=True)

    class Meta:
        filtering = {
            "std_ctab": ALL_WITH_RELATIONS,
            "ctab": ALL,
            "multiple_batch_id": ALL_WITH_RELATIONS,
            "project": ALL_WITH_RELATIONS,
            "with_substructure": ALL_WITH_RELATIONS,
            "similar_to": ALL_WITH_RELATIONS,
            "flexmatch": ALL_WITH_RELATIONS,
            "created": ['gte', 'lte'],
            "created_by": ALL_WITH_RELATIONS,
            "excludes": ALL_WITH_RELATIONS,
        }
        always_return_data = True
        prefix = "related_molregno"
        fieldnames = [('chembl_id', 'chemblId'),
                      ('pref_name', 'preferredCompoundName'),
                      ('max_phase', 'knownDrug'),
                      ('compoundproperties.med_chem_friendly',
                       'medChemFriendly'),
                      ('compoundproperties.ro3_pass', 'passesRuleOfThree'),
                      ('compoundproperties.full_molformula',
                       'molecularFormula'),
                      ('compoundstructures.canonical_smiles', 'smiles'),
                      ('compoundstructures.standard_inchi_key', 'stdInChiKey'),
                      ('compoundproperties.molecular_species', 'species'),
                      ('compoundproperties.num_ro5_violations',
                       'numRo5Violations'),
                      ('compoundproperties.rtb', 'rotatableBonds'),
                      ('compoundproperties.mw_freebase', 'molecularWeight'),
                      ('compoundproperties.alogp', 'alogp'),
                      ('compoundproperties.acd_logp', 'acdLogp'),
                      ('compoundproperties.acd_logd', 'acdLogd'),
                      ('compoundproperties.acd_most_apka', 'acdAcidicPka'),
                      ('compoundproperties.acd_most_bpka', 'acdBasicPka'),
                      ('compoundproperties.full_molformula', 'fullMolformula'),
                      ('project.name', 'Project'),
                        ('multiple_batch_id', 'Upload ID')

                      ]
                         
        fields_to_keep = {
                          'Project': 'Project',
                            'chemblId': 'UOx ID',
                          'id': 'Batch ID',
                          'canonical_smiles': 'SMILES',
                          'created_by': 'Added By',
                          'standard_inchi': 'Std InChi',
                          'molecularWeight': 'Mol Weight',
                          'molecularFormula': 'Mol Formula',
                          'custom_fields': 'custom_fields',
                          'multiple_batch_id' : 'Upload ID' }
        ordrered_ftk = OrderedDict([('chemblId', 'UOx ID'),
                                    ('canonical_smiles', 'SMILES'),
                                    ('knownDrug', 'Known Drug'),
                                    ('medChemFriendly', 'MedChem Friendly'),
                                    ('standard_inchi', 'Std InChi'),
                                    ('molecularWeight', 'Mol Weight'),
                                    ('acdLogp', 'alogp')])

        queryset = CBHCompoundBatch.objects.all()
        resource_name = 'cbh_compound_batches'
        authorization = ProjectAuthorization()
        include_resource_uri = False
        always_return_data = True
        serializer = CBHCompoundBatchSerializer()
        allowed_methods = ['get', 'post', 'put', 'patch']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        paginator_class = Paginator
        # elasticsearch config items
        es_index_name = "chemreg_chemical_index"

    def apply_filters(self, request, applicable_filters):
        """
        An ORM-specific implementation of ``apply_filters``.
        The default simply applies the ``applicable_filters`` as ``**kwargs``,
        but should make it possible to do more advanced things.
        """
        pr = request.GET.get("project__project_key__in", None)
        if pr == "":
            applicable_filters.pop("project__project_key__in")
        #copy reques.GET so we can pass things in it
        request.GET = request.GET.copy()
        smiles = ""
        ws = request.GET.get("with_substructure", None)
        st = request.GET.get("similar_to", None)
        fm = request.GET.get("flexmatch", None)
        # this is the similarity index for fingerprint-like searching
        fp = request.GET.get("fpValue", None)
        cms = None
        pids = self._meta.authorization.project_ids(request)
        # get filters which indicate blanks should be ignored
        if ws:
            smiles = self.convert_mol_string(ws)
            cms = CompoundMols.objects.with_substructure(smiles)
        elif st:
            smiles = self.convert_mol_string(st)
            if fp == None:
                cms = CompoundMols.objects.similar_to(smiles, 90)
            else:
                cms = CompoundMols.objects.similar_to(smiles, fp)
        elif fm:
            smiles = self.convert_mol_string(fm)
            cms = CompoundMols.objects.flexmatch(smiles)
        request.GET["substructure_smarts"] = smiles
        # else:
        #    cms = CompoundMols.objects.all()
        # To be generalised
        cust = request.GET.get("search_custom_fields__kv_any", None)
        # initialise this with project ids, sets up the correct AND
        # initialisation with custom fields
        '''
        cust_queries = Q(project_id__in=set(pids))
        if cust:
            # loop through the custom fields
            # split on pipe (|)
            # put items which are from the same custom field into an OR Q query
            # https://docs.djangoproject.com/en/1.7/topics/db/queries/#complex-lookups-with-q-objects
            cfields = json.loads(cust)
            grouped_fields = {}
            for cfield in cfields:
                cfield_parts = cfield.split("|")
                if grouped_fields.has_key(cfield_parts[0]):
                    grouped_fields[cfield_parts[0]].append(cfield)
                else:
                    grouped_fields[cfield_parts[0]] = [cfield]

            #grouped_fields = json.dumps(grouped_fields)
            for key, val in grouped_fields.iteritems():
                print val
                field_specific_queries = [
                    Q(custom_fields__kv_single=value) for value in val]
                # initialise with the first object
                inner_queries = field_specific_queries.pop()
                for item in field_specific_queries:
                    # OR the subqueries from the same custom field column
                    inner_queries |= item
                # AND the sets of custom field queries
                cust_queries &= inner_queries
            #applicable_filters["custom_fields__kv_any"] = cust
        '''
        if cms != None:
            # run the sql for pulling in new compounds into compound_mols
            indexed = CBHCompoundBatch.objects.index_new_compounds()
            applicable_filters["related_molregno_id__in"] = cms.values_list(
                "molecule_id", flat=True)
        if request.GET.get("related_molregno__chembl__chembl_id__in", None):
            applicable_filters["related_molregno__chembl__chembl_id__in"] = request.GET.get(
                "related_molregno__chembl__chembl_id__in").split(",")

        dateend = applicable_filters.get("created__lte", None)
        if dateend:
            applicable_filters["created__lte"] += " 23:59:59"
        dataset = self.get_object_list(request).filter(
            **applicable_filters)     #.filter(cust_queries) done in elasticsearch instead
        func_group = request.GET.get("functional_group", None)

        if func_group:
            funccms = CompoundMols.objects.with_substructure(func_group)
            dataset = dataset.filter(
                related_molregno_id__in=funccms.values_list("molecule_id", flat=True))
        return dataset.order_by("-id")

    def get_chembl_ids(self, request, **kwargs):
        '''Get a single list of pinned fields for the project previously listed all custom fields in the DB but this was unwieldy'''
        bundle = self.build_bundle(request=request)
        pids = self._meta.authorization.project_ids(request)
        filters = {"project__id__in": pids}
        prefix = request.GET.get(
            "chembl_id__chembl_id__startswith", None).upper()
        desired_format = self.determine_format(request)
        if(prefix):
            filters["chembl_id__chembl_id__startswith"] = prefix
            uox_ids = list(MoleculeDictionary.objects.filter(
                **filters).values_list("chembl_id", flat=True)[0:20])
            bundle.data = [{"value": uox, "label": uox} for uox in uox_ids]
            serialized = json.dumps(bundle.data)
        else:
            serialized = "[]"

        rc = HttpResponse(
            content=serialized, content_type=build_content_type(desired_format), )
        return rc

    def get_elasticsearch_ids(self, request, **kwargs):
        '''Fetch the UOx ids from elasticsearch'''
        bundle = self.build_bundle(request=request)
        pids = self._meta.authorization.project_ids(request)
        filters = {"project__id__in": pids}
        prefix = request.GET.get(
            "chembl_id__chembl_id__startswith", None).upper()
        desired_format = self.determine_format(request)
        if(prefix):
            #filters["chembl_id__chembl_id__startswith"] = prefix
            #uox_ids = list(MoleculeDictionary.objects.filter(**filters).values_list("chembl_id", flat=True)[0:20])
            uox_ids = list(
                elasticsearch_client.get_autocomplete(pids, prefix, 'chemblId'))
            blinded_ids = list(elasticsearch_client.get_autocomplete(pids, prefix, 'blindedBatchId'))
            bundle.data = [{"value": uox, "label": uox} for uox in uox_ids + blinded_ids]
            serialized = json.dumps(bundle.data)
        else:
            serialized = "[]"

        rc = HttpResponse(
            content=serialized, content_type=build_content_type(desired_format), )
        return rc

    def get_elasticsearch_autocomplete(self, request, **kwargs):
        '''Fetch any autocomplete data from elasticsearch'''
        bundle = self.build_bundle(request=request)
        pids = self._meta.authorization.project_ids(request)
        prefix = request.GET.get("custom__field__startswith", None)
        custom_field = request.GET.get("custom_field", None)
        desired_format = self.determine_format(request)
        # send these project ids to the elasticsearch query?
        #filters["search_custom_fields__kv_any"] = prefix
        uox_ids = list(elasticsearch_client.get_autocomplete(pids,
                                                             prefix,
                                                             'custom_field_list.aggregation',
                                                             custom_fields=True,
                                                             single_field=custom_field))
        #bundle.data = ["%s|%s" % (uox, uox) for uox in uox_ids]
        bundle.data = [{"value": uox, "label": self.labelify_aggregate(
            uox, single_field=custom_field)} for uox in uox_ids]
        serialized = json.dumps(bundle.data)
        rc = HttpResponse(
            content=serialized, content_type=build_content_type(desired_format), )
        return rc

    def labelify_aggregate(self, agg, single_field=None):
        # remove the pipe and add square brackets to the custom field type
        splits = agg.split("|")
        label = agg
        # if this is for a single field or custom field, we don't need to show
        # the facet name
        if(len(splits) > 1 and single_field is None):
            label = '%s: %s' % (splits[0], splits[1])
        elif(len(splits) > 1):
            label = splits[1]
        return label

    def reindex_elasticsearch(self, request, **kwargs):
        desired_format = self.determine_format(request)
        batches = self.get_object_list(request)
        # we only want to store certain fields in the search index
        from django.core.paginator import Paginator
        paginator = Paginator(batches, 1000) # chunks of 1000

        for page in range(1, paginator.num_pages +1):
            bs = paginator.page(page).object_list
            batch_dicts = self.batches_to_es_ready(
                bs, request, non_chem_data_only=True)
            # reindex compound data
            index_name = elasticsearch_client.get_main_index_name()
            es_reindex = elasticsearch_client.create_temporary_index(
                batch_dicts, request, index_name)
            if es_reindex.get("errors"):
                print "ERRORS"
                print json.dumps(es_reindex)
            # here you can do what you want with the row
            
            print "done page %d of %d" % (page ,paginator.num_pages)

       
        return HttpResponse(content="test", content_type=build_content_type(desired_format))

    def reindex_compound(self, request, **kwargs):
        # call this when we need to re-index a compound record which has had
        # fields edited
        desired_format = self.determine_format(request)
        deserialized = self.deserialize(
            request, request.body, format=desired_format)
        bundle = self.build_bundle(
            data=dict_strip_unicode_keys(deserialized), request=request)
        id = bundle.data.get("id", None)
        dataset = self.get_object_list(request).filter(id=id)
        batch_dicts = self.batches_to_es_ready(
            dataset, request, non_chem_data_only=True)
        es_reindex = elasticsearch_client.reindex_compound(batch_dicts[0], id)
        return HttpResponse(content=json.dumps({"data": es_reindex}), content_type=build_content_type(desired_format))

    def convert_mol_string(self, strn):
        # commit
        try:
            mol = Chem.MolFromMolBlock(strn)
            smiles = Chem.MolToSmiles(mol)
        except:
            smiles = strn
        self.substructure_smarts = smiles
        return smiles

    def match_list_to_moleculedictionaries(self, batch, project, structure_type="MOL"):
        if structure_type == "MOL":
            structure_key = batch.standard_inchi_key
        else:
            raise NotImplemented
        proj_data = MoleculeDictionary.objects.by_project_and_natural_key(structure_type,
                                                                          structure_key,
                                                                          project.pk)
        same_project = proj_data.values("molregno", "chembl", "created_by")

        pub = MoleculeDictionary.objects.by_natural_key_public_except_project(
            structure_type,
            structure_key,
            project.pk).order_by("-insert_date")
        different_project_but_public = pub.values(
            "molregno", "chembl", "created_by")
        all_items = list(same_project) + list(different_project_but_public)
        linkedproject = 0
        linkedpublic = 0
        new = 0
        for item in all_items:
            item["tobelinked"] = False
        if same_project.count() > 0:
            # Increment the forced registration number compared to what is
            # already in the database as this can then be used to force the
            # registration of the molecule
            forced_reg_no = proj_data.aggregate(
                Max('forced_reg_index'))["forced_reg_index__max"] + 1
            linkedproject += 1
            batch.warnings["forced_reg_no"] = forced_reg_no
        elif different_project_but_public.count() > 0:
            linkedpublic += 1
        else:
            new += 1
        batch.warnings["linkedpublic"] = linkedpublic
        batch.warnings["linkedproject"] = linkedproject
        if len(all_items) > 0:
            all_items[0]["tobelinked"] = True
        batch.warnings["linkable_molecules"] = all_items

    def post_validate(self, request, **kwargs):
        """Runs the validation for a single or small set of molecules"""
        #self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
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
        #updated_bundle = self.obj_build(bundle, dict_strip_unicode_keys(deserialized))
        try:
            m = Chem.MolFromMolBlock(bundle.data["ctab"])
            if not m:
                raise Exception("valancy_or_other_error")
        except:
            raise Exception("valancy_or_other_error")
        obj = CBHCompoundBatch.objects.from_rd_mol(
            m, orig_ctab=bundle.data["ctab"], )
        bundle.obj = obj
        bundle.obj.validate()
        self.match_list_to_moleculedictionaries(
            bundle.obj, bundle.data["project"])
        dictdata = bundle.obj.__dict__
        dictdata.pop("_state")

        updated_bundle = self.build_bundle(obj=bundle.obj, data=dictdata)
        return self.create_response(request, updated_bundle, response_class=http.HttpAccepted)

    def save_related(self, bundle):
        
        bundle.obj.generate_structure_and_dictionary()

    def alter_deserialized_list_data(self, request, deserialized):
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
        bundle = super(CBHCompoundBatchResource, self).full_hydrate(bundle)
        if not bundle.obj.id:
            bundle.obj.created_by = bundle.request.user.username
            bundle.obj.created_by_id = bundle.request.user.id
            # bundle.obj.validate()
            # self.match_list_to_moleculedictionaries(bundle.obj,bundle.data["project"] )
        return bundle

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/delete_index/?$" % self._meta.resource_name,
                self.wrap_view('delete_index'), name="delete_index"),
            url(r"^(?P<resource_name>%s)/update_temp_batches/?$" % self._meta.resource_name,
                self.wrap_view('update_temp_batches'), name="update_temp_batches"),
            url(r"^(?P<resource_name>%s)/get_part_processed_multiple_batch/?$" % self._meta.resource_name,
                self.wrap_view('get_part_processed_multiple_batch'), name="api_get_part_processed_multiple_batch"),
            url(r"^(?P<resource_name>%s)/get_list_elasticsearch/?$" % self._meta.resource_name,
                self.wrap_view('get_list_elasticsearch'), name="api_get_list_elasticsearch"),
            url(r"^(?P<resource_name>%s)/get_chembl_ids/?$" % self._meta.resource_name,
                self.wrap_view('get_chembl_ids'), name="api_get_chembl_ids"),
            url(r"^(?P<resource_name>%s)/get_elasticsearch_ids/?$" % self._meta.resource_name,
                self.wrap_view('get_elasticsearch_ids'), name="api_get_elasticsearch_ids"),
            url(r"^(?P<resource_name>%s)/reindex_elasticsearch/?$" % self._meta.resource_name,
                self.wrap_view('reindex_elasticsearch'), name="api_compounds_reindex_elasticsearch"),
            url(r"^(?P<resource_name>%s)/reindex_compound/?$" % self._meta.resource_name,
                self.wrap_view('reindex_compound'), name="api_reindex_compound"),
            url(r"^(?P<resource_name>%s)/get_elasticsearch_autocomplete/?$" % self._meta.resource_name,
                self.wrap_view('get_elasticsearch_autocomplete'), name="api_get_elasticsearch_autocomplete"),
            url(r"^(?P<resource_name>%s)/validate/?$" % self._meta.resource_name,
                self.wrap_view('post_validate'), name="api_validate_compound_batch"),
            url(r"^(?P<resource_name>%s)/validate_list/?$" % self._meta.resource_name,
                self.wrap_view('post_validate_list'), name="api_validate_compound_list"),
            url(r"^(?P<resource_name>%s)/validate_drawn/?$" % self._meta.resource_name,
                self.wrap_view('post_validate_drawn'), name="api_validate_compound_drawn"),
            url(r"^(?P<resource_name>%s)/multi_batch_save/?$" % self._meta.resource_name,
                self.wrap_view('multi_batch_save'), name="multi_batch_save"),
            url(r"^(?P<resource_name>%s)/multi_batch_custom_fields/?$" % self._meta.resource_name,
                self.wrap_view('multi_batch_custom_fields'), name="multi_batch_custom_fields"),
            url(r"^(?P<resource_name>%s)/validate_files/?$" % self._meta.resource_name,
                self.wrap_view('post_validate_files'), name="api_compound_validate_files"),
            url(r"^(?P<resource_name>%s)/export_file/?$" % self._meta.resource_name,
                self.wrap_view('export_file'), name="api_compound_export_file"), 
            url(r"^(?P<resource_name>%s)/get_saved_searches/?$" % self._meta.resource_name,
                self.wrap_view('get_saved_search'), name="api_compound_get_saved_searches"), 
            url(r"^(?P<resource_name>%s)/post_saved_searches/?$" % self._meta.resource_name,
                self.wrap_view('post_saved_search'), name="api_compound_post_saved_searches"), ]


    def get_saved_searches(self, request, **kwargs):
        pass

    def post_saved_searches(self, request, **kwargs):
        pass

    def patch_dict(self, dictdata, headers):
        for header in headers:
            json_patches = copy.copy(header.get("operations", False))
            if json_patches:
                apply_json_patch(dictdata, json_patches)

    def multi_batch_save(self, request, **kwargs):
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
        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        project = Serializer().serialize(mb)
        limit = 500
        offset = 0
        batches = []
        hasMoreData = True
        while hasMoreData:

            bundles = self.get_cached_temporary_batch_data(
                mb.id,  {"limit": limit, "offset": offset, "sorts": '[{"id": {"order": "desc"}}]'}, request)
            #allowing setting of headers to be fale during saving of drawn molecule
            if bundle.data["headers"]:
                for d in bundles["objects"]:
                    self.patch_dict(d, copy.deepcopy(bundle.data["headers"]))
            set_of_batches = self.get_cached_temporary_batches(
                bundles, request,  bundledata=bundle.data)
            batches.extend(set_of_batches["objects"])
            offset += limit
            if len(set_of_batches["objects"]) == 0:
                hasMoreData = None
        mb.saved = True
        mb.created_by = request.user.username
        mb.save()
        bundle.data["saved"] = 0
        bundle.data["ignored"] = 0
        to_be_saved = []
        for batch in batches:
            if(batch.obj.properties.get("action", "") == "New Batch"):
                batch.obj.created_by = request.user.username
                batch.obj.created_by_id = request.user.id
                batch.obj.id = None
                batch.obj.multiple_batch_id = id
                batch.obj.generate_structure_and_dictionary()

                batch.multi_batch_id = id
                bundle.data["saved"] += 1
                to_be_saved.append(batch.obj)
            index_batches_in_new_index(to_be_saved)
        if(mb.uploaded_file):
            if(mb.uploaded_file.extension == ".sdf"):
                newreq = copy.copy(request)
                self.alter_batch_data_after_save(
                    to_be_saved, 
                    mb.uploaded_file.file,
                    newreq,
                    mb
                )
        elasticsearch_client.delete_index(
            elasticsearch_client.get_temp_index_name(request, mb.id))
        batch_dicts = self.batches_to_es_ready(to_be_saved, request)
        index_name = elasticsearch_client.get_main_index_name()
        elasticsearch_client.create_temporary_index(
            batch_dicts, request, index_name)
        self.after_save_and_index_hook(request, id, mb.project.project_key)
        # this needs to be the main elasticsearch compound index
        # and should update any existing records in there? That might be in
        # another method
        return self.create_response(request, bundle, response_class=http.HttpCreated)

    def after_save_and_index_hook(self, request, multi_batch_id, project_key):

        pass

    def alter_batch_data_after_save(self, batch_list, python_file_object,request, multi_batch):
        pass


    def delete_index(self, request, **kwargs):
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
        mb = CBHCompoundMultipleBatch.objects.get(pk=id)
        elasticsearch_client.delete_index(
            elasticsearch_client.get_temp_index_name(request, mb.id))
        return self.create_response(request, bundle, response_class=http.HttpAccepted)

    def convert_custom_field(self, uncurated_value, field_schema):
        curated_value = uncurated_value
        if field_schema.get("format", "") == "yyyy-mm-dd":
            if uncurated_value:

                curated_value = dateutil.parser.parse(
                    uncurated_value).strftime("%Y-%m-%d")
        elif (field_schema.get("field_type", "") == PinnedCustomField.UISELECTTAGS):

            curated_value = json.dumps(uncurated_value.split(","))
        elif (field_schema.get("field_type", "") == PinnedCustomField.INTEGER):
            curated_value = int(uncurated_value)
        elif (field_schema.get("field_type", "") in [PinnedCustomField.NUMBER, PinnedCustomField.PERCENTAGE]):
            curated_value = float(uncurated_value)
        return curated_value

    def update_temp_batches(self, request, **kwargs):
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
        multi_batch_id = bundle.data["multiplebatch"]
        es_serializer = CBHCompoundBatchElasticSearchSerializer()
        es_ready_updates = [es_serializer.to_es_ready_data(dictdata,
                                                           options={"underscorize": True}) for dictdata in bundle.data["objects"]]
        index_name = elasticsearch_client.get_temp_index_name(
            request, multi_batch_id)
        elasticsearch_client.create_temporary_index(
            es_ready_updates, request, index_name)
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

    def validate_multi_batch(self, multi_batch, bundle, request, batches):
        batches_not_errors = [batch for batch in batches if batch and not batch.warnings.get(
            "parseerror", None) and not batch.warnings.get("smilesParseError", None)]

        if (len(batches) == 0):
            raise BadRequest("no_data")
        for b in batches_not_errors:
            b.properties["action"] = "New Batch"
        batches_with_structures = [
            batch for batch in batches_not_errors if not batch.blinded_batch_id]
        blinded_data = [
            batch for batch in batches_not_errors if batch.blinded_batch_id]
        sdfstrings = [batch.ctab for batch in batches_with_structures]
        sdf = "\n".join(sdfstrings)

        filename = "/tmp/" + shortuuid.ShortUUID().random()
        text_file = open(filename, "w")
        text_file.write(sdf)
        text_file.close()
        from subprocess import PIPE, Popen
        p = Popen([settings.INCHI_BINARIES_LOCATION['1.02'],
                   "-STDIO",  filename], stdout=PIPE, stderr=PIPE)
        a = p.communicate()
        inchis = {}

        #PB - there is an assumption here that everything that has a structure will generate an inChi without issue. This is not the case.
        #Where a molecule does not generate an inchi, there will be a key error looking up the inchi in inchiparts, as anything that cannot 
        #generate an inchi will be missing from inchiparts, i.e. 50 structures with 1 error will have 49 entries in inchiparts, and this 
        #will in turn bin the whole file - not great when we can handle erroring structures elsewhere

        error_locs = []

        #a[0] holds the generated inchis. a[1] holds all of the error and warning information (if any)
        errorparts = a[1].split("\nError")
        if(len(errorparts) > 1):
            for i, errorp in enumerate(errorparts):
                #split on 'structure #', then get the number given
                if(i > 0):
                    splits = errorp.split('structure #')

                    error_loc = splits[1].split('.')[0]
                    #convert to number, put this number in an errors list
                    error_locs.append(error_loc)

        err_batches = []
        #for the errors found, remove from non-error lists and flag as erroring
        for error_no in error_locs:
            error_no_int = int(float(error_no)) - 1

            #find structures at the position indicated - 1 (for 0-indexed list)
            err_batch = batches_with_structures[error_no_int]
            err_batches.append(err_batch)
 
        #we can't remove these while looping through err_locs as it messes up the list order and gives arrayindex exceptions
        for err_batch in err_batches:

            #remove from batches_with_structures and batches_not_errors
            batches_with_structures.remove(err_batch)
            batches_not_errors.remove(err_batch)

            #flag this batch as erroring due to inability to generate anything for the standard_inchi_key field
            batches_index = batches.index(err_batch)
            batches[batches_index].warnings["inchicreationerror"] = "true"
            batches[batches_index].properties["action"] = "Ignore"


        inchiparts = a[0].split("\nStructure:")

        for i, inch in enumerate(inchiparts):
            parts = inch.split("\n")
            if len(parts) == 1:
                continue
            ints = [s for s in parts[0].split() if s.isdigit()]
            part = "".join(ints)
            inchis[part] = parts[1]
        if not bundle.data.get("fileerrors"):
            bundle.data["fileerrors"] = []
        new_uploaded_data = []
        already_found = set([])
        duplicates = set([])
        for i, batch in enumerate(batches_with_structures):
            if (str(i+1) in error_locs):
                batch.standard_inchi = None
            else: 
                batch.standard_inchi = inchis[str(i+1)]
            batch.validate(temp_props=False)
            if batch.standard_inchi_key in already_found:
                # setting this in case we change it later
                duplicates.add(batch.standard_inchi_key)
            else:
                already_found.add(batch.standard_inchi_key)

            new_uploaded_data.append(batch)
        already_in_db = MoleculeDictionary.objects.filter(project=bundle.data[
                                                          "project"], structure_type="MOL", structure_key__in=already_found).values_list("structure_key", flat=True)
        already_in_db = set(already_in_db)

        bundle.data["new"] = 0
        new_data = set([])
        duplicate_overlaps = set([])
        duplicate_new = set([])
        for batch in batches_with_structures:
            if batch.standard_inchi_key in duplicates:
                batch.warnings["duplicate"] = True
            if batch.standard_inchi_key in already_in_db:
                batch.warnings["overlap"] = True
                if batch.standard_inchi_key in duplicates:
                    batch.warnings["duplicate"] = True
                    duplicate_overlaps.add(batch.standard_inchi_key)
            else:
                batch.warnings["new"] = True

                new_data.add(batch.standard_inchi_key)
                if batch.standard_inchi_key in duplicates:
                    batch.warnings["duplicate"] = True
                    duplicate_new.add(batch.standard_inchi_key)

        for batch in batches_with_structures:
            if batch.warnings.get("nostructure") == True:
                del batch.warnings["nostructure"]
        for batch in blinded_data:
            batch.warnings["nostructure"] = True
        bundle.data["batchstats"] = {}
        bundle.data["batchstats"]["withstructure"] = len(
            batches_with_structures)
        bundle.data["batchstats"]["parseerrors"] = len(batches) - len(batches_not_errors) + len(
            [b for b in batches_not_errors if b.warnings.get("parseerror", False) == "true"])
        bundle.data["batchstats"]["withoutstructure"] = len(blinded_data)
        bundle.data["batchstats"]["total"] = len(batches)
        bundle.data["compoundstats"] = {}
        bundle.data["compoundstats"]["total"] = len(
            already_in_db) + len(new_data)
        bundle.data["compoundstats"]["overlaps"] = len(already_in_db)
        bundle.data["compoundstats"]["new"] = len(new_data)
        bundle.data["compoundstats"][
            "duplicateoverlaps"] = len(duplicate_overlaps)
        bundle.data["compoundstats"]["duplicatenew"] = len(duplicate_new)
        bundle.data["multiplebatch"] = multi_batch.pk

        
        fifty_batches_for_first_page = self.set_cached_temporary_batches(
            batches, multi_batch.id, request)
        multi_batch.uploaded_data = bundle.data
        multi_batch.save()
        #bundle.data["objects"] = fifty_batches_for_first_page
        index_name = elasticsearch_client.get_temp_index_name(
            request, multi_batch.id)
        elasticsearch_client.get_action_totals(index_name, bundle.data)
        return self.create_response(request, bundle, response_class=http.HttpAccepted)

    def get_part_processed_multiple_batch(self, request, **kwargs):
        """
        Get the part processed data from elasticsearch and the stats about the
        multiple batch
        """
        # TODO: Uncached for now. Invalidation that works for everyone may be
        #       impossible.
        bundle = self.build_bundle(request=request)

        # self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        if(kwargs.get("multi_batch", None)):
            mb = kwargs.get("multi_batch")
            id = mb.id
        else:
            id = request.GET.get("current_batch")
            mb = CBHCompoundMultipleBatch.objects.get(pk=id)

        to_be_serialized = mb.uploaded_data
        to_be_serialized = self.get_cached_temporary_batch_data(
            id, request.GET, request, bundledata=to_be_serialized)
        index_name = elasticsearch_client.get_temp_index_name(request, id)
        elasticsearch_client.get_action_totals(index_name, to_be_serialized)
        to_be_serialized = self.alter_list_data_to_serialize(request, to_be_serialized, permanent_data=False)
        return self.create_response(request, to_be_serialized)

    def post_validate_list(self, request, **kwargs):
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
        multiple_batch = CBHCompoundMultipleBatch.objects.create(
            project=bundle.data["project"])
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
            batches.append(b)

        bundle.data["current_batch"] = multiple_batch.pk
        bundle.data["headers"] = []
        return self.validate_multi_batch(multiple_batch, bundle, request, batches)

    def post_validate_drawn(self, request, **kwargs):
        """
        Turn a drawn molfile via chemdoodle into a processable rdkit mol to be used by the rest of the architecture
        """
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
        #now get the SDfile from the request
        ctab = bundle.data.get('sketch', '')
        # first assume ctab
        # allmols = [(obj, Chem.MolFromSmiles(str(obj))) for obj in objects]
        # we need to make allmols a single-entry list in order to use the rst of the architecture

        mol = Chem.MolFromMolBlock(ctab)
        #for each of the supplied custom fields, use set prop to replicate how this would be if it was an sd file.
        prj_data_fields = bundle.data.get('custom_fields',{})
        sup_data_fields = bundle.data.get('supplementary_fields',[])
        if(mol):
            for k1, v1 in prj_data_fields.iteritems():
                if(isinstance(v1, list)):
                    v1str = ""
                    for item in v1:
                        v1str = v1str + ", " + item.encode('ascii', 'ignore')
                    mol.SetProp(k1.encode('ascii', 'ignore'), v1str)
                #not a unicode string? Don't try and ascii encode the result, 
                #just add it as a value converted to a string (which SetProp expects)
                elif(isinstance(v1, basestring) != True):
                    mol.SetProp(k1.encode('ascii', 'ignore'), str(v1))
                else:
                    mol.SetProp(k1.encode('ascii', 'ignore'), v1.encode('ascii', 'ignore'))

            for field in sup_data_fields:
                mol.SetProp(field['name'].encode('ascii', 'ignore'), field['value'].encode('ascii', 'ignore'))
        
        allmols = [ mol ]
        errors = []
        batches = []
        multiple_batch = CBHCompoundMultipleBatch.objects.create(
            project=bundle.data["project"])
        try:
            b = CBHCompoundBatch.objects.from_rd_mol(
                mol, project=bundle.data["project"])
            batches.append(b)
        except Exception, e:
            b = CBHCompoundBatch.objects.blinded(
                            project=bundle.data["project"])
            b.warnings["parseerror"] = "true"
            b.properties["action"] = "Ignore"
            errors.append(
                    {"message": str(e)})
        b.custom_fields = prj_data_fields
        b.uncurated_fields = {}
        #iterate the supplementary fields
        for field in sup_data_fields:
            b.uncurated_fields[field['name'].encode('ascii', 'ignore')] = field['value'].encode('ascii', 'ignore')
        batches.append(b)
        bundle.data["current_batch"] = multiple_batch.pk
        bundle.data["errors"] = errors
        bundle.data["headers"] = False
        return self.validate_multi_batch(multiple_batch, bundle, request, batches)

    def preprocess_sdf_file(self, file_obj, request):
        pass

    def post_validate_files(self, request, **kwargs):

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
        correct_file = FlowFile.objects.get(
            identifier="%s-%s" % (session_key, file_name))
        batches = []
        headers = []
        errors = []
        fielderrors = {}
        fielderrors["stringdate"] = set([])
        fielderrors["number"] = set([])
        fielderrors["integer"] = set([])
        structure_col = bundle.data.get("struccol", "")
        if (".cdx" in correct_file.extension):
            mols = [mol for mol in readfile(
                str(correct_file.extension[1:]), str(correct_file.file.name), )]
            rxn = None
            index = 0
            if correct_file.extension == '.cdxml':
                # Look for a stoichiometry table in the reaction file
                rxn = chemdraw_reaction.parse(str(correct_file.file.name))
                headers = ["%Completion",
                           "%Yield",
                           "Expected Moles",
                           "Product Moles",
                           "Expected Mass",
                           "Product Mass",
                           "MW",
                           "role",
                           "Purity",
                           "Limit Moles",
                           # "Formula",
                           "Equivalents",
                           "Measured Mass"]

            for pybelmol in mols:
                molfile = pybelmol.write("mdl")
                if molfile.strip() and molfile != "*":
                    rd_mol = Chem.MolFromMolBlock(molfile, sanitize=False)
                    '''
                    Pybel can read some hypervalent molecules that RDKit cannot read
                    Therefore currently these molecules are outputted as images and sent back to the front end
                    https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg04466.html
                    '''
                    if rd_mol:
                        smiles = Chem.MolToSmiles(rd_mol)
                        if smiles.strip():
                            try:
                                b = CBHCompoundBatch.objects.from_rd_mol(
                                    rd_mol, orig_ctab=molfile, smiles=smiles, project=bundle.data["project"])
                            except Exception, e:
                                b = None
                                index = index - 1
                            if b:
                                if rxn:
                                    # Here we set the uncurated fields equal to
                                    # the reaction data extracted from Chemdraw
                                    b.uncurated_fields = rxn.get(
                                        pybelmol.title, {})
                                batches.append(b)
                            else:
                                errors.append({"index": index+1, "image": pybelmol.write(
                                    "svg"), "message": "Unable to produce inchi from this molecule"})
                    else:
                        errors.append({"index": index+1, "image": pybelmol.write(
                            "svg"), "message": "Invalid valency or other error parsing this molecule"})
                    index += 1

        else:
            if (correct_file.extension == ".sdf"):
                # read in the file
                self.preprocess_sdf_file(correct_file.file, request)
                suppl = Chem.ForwardSDMolSupplier(correct_file.file)
                mols = [mo for mo in suppl]
                if(len(mols) > 10000):
                    raise BadRequest("file_too_large")
                # read the headers from the first molecule

                headers = get_all_sdf_headers(correct_file.file.name)
                
                uncurated, ctabs = get_uncurated_fields_from_file(correct_file, fielderrors)
                for index, mol in enumerate(mols):
                    if mol is None:
                        b = CBHCompoundBatch.objects.blinded(
                            project=bundle.data["project"])
                        b.warnings["parseerror"] = "true"
                        b.properties["action"] = "Ignore"

                        b.uncurated_fields = uncurated[index]
                        errors.append(
                                {"index": index+1,  "message": "No structure found"})
                    else:
                        try:
                            b = CBHCompoundBatch.objects.from_rd_mol(
                                mol, orig_ctab=ctabs[index], project=bundle.data["project"])
                        except Exception, e:
                            b = CBHCompoundBatch.objects.blinded(
                            project=bundle.data["project"])
                            b.warnings["parseerror"] = "true"
                            b.properties["action"] = "Ignore"
                            errors.append(
                                    {"index": index+1,  "message": str(e)})
                       
                        b.uncurated_fields = uncurated[index]

                    batches.append(b)

            elif(correct_file.extension in (".xls", ".xlsx")):
                # we need to know which column contains structural info - this needs to be defined on the mapping page and passed here
                # read in the specified structural column

                headerswithdata = set([])

                df = None
                try:
                    df = pd.read_excel(correct_file.file)
                except IndexError:
                    raise BadRequest("no_headers")
                if len(df.index) > 2000:
                    raise BadRequest("file_too_large")
                # read the smiles string value out of here, when we know which
                # column it is.
                row_iterator = df.iterrows()
                headers = list(df)
                headers = [h.replace(".", "__") for h in headers]
                df.columns = headers
                # Only automap on the first attempt at mapping the smiles
                    # column
                if not structure_col and not bundle.data.get("headers", None):
                    max_score = 0
                    for header in headers:
                        # fuzzy matching for smiles - this should also
                        # match things like "canonical_smiles"
                        hdr = re.sub('[^0-9a-zA-Z]+', ' ', header)
                        for h in hdr.split(" "):
                            h = h.strip()
                            if h:
                                score = fuzzymatch(
                                    a="smiles", b=h.lower()).ratio()
                                if score > max_score and score > 0.9:
                                    structure_col = header
                                    max_score = score
                                    automapped_structure = True
                
                for index, row in row_iterator:
                    
                    if structure_col:
                        smiles_str = row[structure_col]
                        try:
                            struc = Chem.MolFromSmiles(smiles_str)
                            if struc:
                                Compute2DCoords(struc)
                                b = CBHCompoundBatch.objects.from_rd_mol(
                                    struc, smiles=smiles_str, project=bundle.data["project"], reDraw=True)
                                b.blinded_batch_id = None
                            else:
                                raise Exception("Smiles not processed")

                        except Exception, e:
                            b = CBHCompoundBatch.objects.blinded(
                                project=bundle.data["project"])
                            b.original_smiles = smiles_str
                            b.warnings["parseerror"] = "true"
                            b.properties["action"] = "Ignore"
                    else:
                        b = CBHCompoundBatch.objects.blinded(
                            project=bundle.data["project"])

                    if dict(b.uncurated_fields) == {}:
                            # Only rebuild the uncurated fields if this has not
                            # been done before
                        parse_pandas_record(
                                headers, b, "uncurated_fields", row, fielderrors, headerswithdata)
                   
                    batches.append(b)
                headers = [hdr for hdr in headers if hdr in headerswithdata]
            else:
                raise BadRequest("file_format_error")
        multiple_batch = CBHCompoundMultipleBatch.objects.create(
                project=bundle.data["project"],
                uploaded_file=correct_file
            )
        for b in batches:
            if b:
                b.multiple_batch_id = multiple_batch.pk
                b.created_by = bundle.request.user.username

        bundle.data["fileerrors"] = errors
        bundle.data["automapped"] = 0
        cfr = ChemRegCustomFieldConfigResource()
        jsondata = json.loads(cfr.get_detail(request,
            pk=bundle.data["project"].custom_field_config_id).content)
        schemaform = [field["edit_form"]["form"][0] for field in jsondata["project_data_fields"]]
        if not bundle.data.get("headers", None):
            bundle.data["headers"] = []
            for header in headers:
                copyto = ""
                automapped = False
                operations = []
                if header == structure_col:
                    copyto = "SMILES for chemical structures"
                    if automapped_structure:
                        automapped = True
                else:
                    form = copy.deepcopy(schemaform)
                    copyto = ""
                    max_score = 0
                    for form_item in form:
                        score = fuzzymatch(
                            a=form_item["key"].lower(), b=header.lower()).ratio()
                        if score > max_score and score > 0.9:
                            matched_item = form_item
                            copyto = matched_item["key"]
                            automapped = True
                            if(matched_item.get("field_type", "") == "uiselecttags"):
                                operations.append(
                                    {"op": "split", "path": "/uncurated_fields/" + header})
                                operations.append(
                                    {"op": "move", "path": "/custom_fields/" + matched_item["key"], "from": "/uncurated_fields/" + header})
                            else:
                                operation = {
                                    "op": "move", "path": "/custom_fields/" + matched_item["key"], "from": "/uncurated_fields/" + header}
                                operations.append(operation)
                                if(matched_item.get("format", "") == "date"):
                                    operations.append(
                                        {"op": "convertdate", "path": "/custom_fields/" + matched_item["key"]})
                            #set the max score so less well matched content than this is ignored
                            max_score = score

                bundle.data["headers"].append({
                    "name": header,
                    "automapped": automapped,
                    "copyto": copyto,
                    "operations": operations,
                    "fieldErrors": {
                        "stringdate": header in fielderrors["stringdate"],
                        "integer": header in fielderrors["integer"],
                        "number": header in fielderrors["number"]
                    }
                })

        return self.validate_multi_batch(multiple_batch, bundle, request, batches)

    def alter_list_data_to_serialize(self, request, data, permanent_data=True):
        '''use the request type to determine which fields should be limited for file download,
           add extra fields if needed (eg images) and enumerate the custom fields into the 
           rest of the calculated fields'''
           #Broken functionality to add images to search results
        # if request.GET.get("substructure_smarts", False):
        #     for index, b in enumerate(data["objects"]):
        #         ctab = b.data["properties"][
        #             "substructureMatch"] = request.GET.get("substructure_smarts")
        if(self.determine_format(request) == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' or request.GET.get("format") == "xlsx"  or request.GET.get("format") == "sdf" or self.determine_format(request) == 'chemical/x-mdl-sdfile'):
            ordered_cust_fields = PinnedCustomField.objects.filter(custom_field_config__project__project_key__in=request.GET.get(
                "project__project_key__in", "").split(",")).order_by("custom_field_config__project__id", "position").values("name", "field_type")
            seen = set()
            deduplicated_cfs = [x for x in ordered_cust_fields if x[
                "name"] not in seen and not seen.add(x["name"])]
            df_data = []
            ordered_cust_fields = []
            # projects = set([])
            uncurated_field_names = set()
            for index, b in enumerate(data["objects"]):
                # remove items which are not listed as being kept
                try:
                    olddata = b.data
                except AttributeError:
                    olddata = b
                if self.determine_format(request) == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' or request.GET.get("format") == "xlsx":
                    #Add the byte encoded image but only if this is an excel export
                    try:
                        new_data = {"Structure Image": b.obj.image}
                    except AttributeError:
                        new_data = {"Structure Image": b.get("image", "")}
                else:
                    new_data = {}
                # projects.add(b.obj.project_id)
                for k, v in olddata.iteritems():
                    for name, display_name in self.Meta.fields_to_keep.iteritems():
                        if k == name:
                            new_data[display_name] = v

                # we need sd format exported results to retain stereochemistry
                # - use mol instaed of smiles
                if(request.GET.get("format") == "sdf"):
                    new_data['ctab'] = olddata['ctab']
                elif (self.determine_format(request) == 'chemical/x-mdl-sdfile' ):
                    new_data['ctab'] = olddata['ctab']
                # dummy
                # not every row has a value for every custom field

                if permanent_data:
                    for item in deduplicated_cfs:
                        cf_value = olddata["custom_fields"].get(item["name"], "")
                        if item["field_type"] == PinnedCustomField.UISELECTTAGS:
                            if isinstance(cf_value, basestring):
                                try:
                                    cf_value = json.loads(cf_value).join(",")
                                except:
                                    pass
                            elif isinstance(cf_value, list):
                                cf_value = cf_value.join(",")
                        new_data[item["name"]] = cf_value

                # now remove custom_fields

                del(new_data['custom_fields'])
                for field, value in olddata['uncurated_fields'].iteritems():
                    new_data[field] = value
                    uncurated_field_names.add(field)
                # #now remove custom_fields
                # del(new_data['uncurated_fields'])
                try:
                    b.data = new_data
                except AttributeError:
                    pass

                df_data.append(new_data)
            df = pd.DataFrame(df_data)
            data['export'] = df.to_json()
            data["headers"] = {
                "cbh": [display_name for name, display_name in self.Meta.fields_to_keep.iteritems()],
                "custom_fields": [cf["name"] for cf in deduplicated_cfs],
                "uncurated_fields": sorted(list(uncurated_field_names))
            }

        return data

    def create_response(self, request, data, response_class=HttpResponse, **response_kwargs):
        """
        Extracts the common "which-format/serialize/return-response" cycle.
        Mostly a useful shortcut/hook.
        """
        
        desired_format = self.determine_format(request)
        serialized = self.serialize(
            request, data, desired_format, options=self.Meta.ordrered_ftk)
        rc = response_class(content=serialized, content_type=build_content_type(
            desired_format), **response_kwargs)
        if(desired_format == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'):
            rc['Content-Disposition'] = 'attachment; filename=export_from_chembiohub_chemreg%d.xlsx' % int(time.time())
        elif(desired_format == 'chemical/x-mdl-sdfile'):
            rc['Content-Disposition'] = 'attachment; filename=export_from_chembiohub_chemreg%d.sdf' % int(time.time())
        return rc

    def dehydrate(self, bundle):
        # try:
        data = bundle.obj.related_molregno
        user = None

        for names in self.Meta.fieldnames:
            if not bundle.obj.blinded_batch_id:
                try:
                    bundle.data[names[1]] = deepgetattr(data, names[0], None)
                except(AttributeError):
                    bundle.data[names[1]] = ""
            else:
                if names[1] == "chemblId":
                    bundle.data[names[1]] = bundle.obj.blinded_batch_id
                else:
                    bundle.data[names[1]] = ""
        if bundle.obj.created_by:
          #user = User.objects.get(username=bundle.obj.created_by)
            User = get_user_model()
            try:
                user = User.objects.get(username=bundle.obj.created_by)
            except ObjectDoesNotExist:
                try:
                    user = User.objects.get(first_name=bundle.obj.created_by.split[" "][
                                            0], last_name=bundle.obj.created_by.split[" "][1])
                except:
                    user = None
        mynames = [ "uncurated_fields", "warnings", "properties", "custom_fields",]
        for name in mynames:
            bundle.data[name] = json.loads(bundle.data[name])
            
        #bundle.data["created_by"] = user.__dict__
        if user != None:
            if user.first_name:
                bundle.data["created_by"] = "%s %s" % (
                    user.first_name, user.last_name)
            else:
                bundle.data["created_by"] = user.username
        else:
            bundle.data["created_by"] = bundle.obj.created_by
        bundle.data["timestamp"] = str(bundle.data["created"])[0:10]
        # except:
        #    pass
        return bundle

    def batches_to_es_ready(self, batches, request, non_chem_data_only=None):
        batch_dicts = []
        index = 1
        es_serializer = CBHCompoundBatchElasticSearchSerializer()
        for batch in batches:
            if batch:
                if(not batch.id):
                    batch.id = index
                bun = self.build_bundle(obj=batch, request=request)
                bun = self.full_dehydrate(bun, for_list=True)
                if bun.data["project"] == '':
                    bun.data[
                        "project"] = '/%s/cbh_projects/%d' % (settings.WEBSERVICES_NAME, bun.obj.project_id)
                if non_chem_data_only:
                    ready = es_serializer.to_es_ready_non_chemical_data(
                        bun.data, options={"underscorize": True})
                else:
                    ready = es_serializer.to_es_ready_data(
                        bun.data, options={"underscorize": True})

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

    def set_cached_temporary_batches(self, batches, multi_batch_id, request):
        es_serializer = CBHCompoundBatchElasticSearchSerializer()
        batch_dicts = self.batches_to_es_ready(batches, request)
        index_name = elasticsearch_client.get_temp_index_name(
            request, multi_batch_id)
        elasticsearch_client.create_temporary_index(
            batch_dicts, request, index_name)
        # Now get rid of my ES preparation again

    def get_cached_temporary_batches(self, bundles, request, bundledata={}):
        fakeobj = CBHCompoundBatch(project=bundledata["project"], id=10)
        bun = self.build_bundle(obj=fakeobj)
        bun = self.full_dehydrate(bun)
        proj = bun.data["project"]
        data = []
        for datum in bundles["objects"]:
            datum = self.build_bundle(data=datum, request=request)
            datum.data["project"] = proj
            datum = self.full_hydrate(datum)
            data.append(datum)
        bundles["objects"] = data
        return bundles

    def get_list_elasticsearch(self, request, **kwargs):
        """
        Search for compound batches or inventory items
        The search GET request is interpretted into a request to elasticsearch.
        If a request for a structure search is made then the search is executed first using PostgreSQL 
        via the RDKit Cartridge. Otherwise the search goes direct to elasticsearch
        """
        # TODO: Uncached for now. Invalidation that works for everyone may be
        #       impossible.
        base_bundle = self.build_bundle(request=request)
        must_list = []
        # If there is a substructure query then we use the non-elasticsearch
        # implementation to pull back the required fields


        ws = request.GET.get("with_substructure", None)
        st = request.GET.get("similar_to", None)
        fm = request.GET.get("flexmatch", None)
        # this is the similarity index for fingerprint-like searching
        fp = request.GET.get("fpValue", None)
        fg = request.GET.get("functional_group", None)


        query = None
        if ws or st or fm or fp or fg:
            objects = self.obj_get_list(
                bundle=base_bundle, **self.remove_api_resource_names(kwargs))
            ids = objects.values_list("id", flat=True)[0:100000]
            query = {
                "ids": {
                "type": "batches",
                "values": [str(i) for i in ids]
                }
            }

        get_data = request.GET
        blanks_filter = json.loads(get_data.get("showBlanks", "[]"))
        blanks_filter = [
            get_key_from_field_name(field) for field in blanks_filter]
        nonblanks_filter = json.loads(get_data.get("showNonBlanks", "[]"))
        nonblanks_filter = [
            get_key_from_field_name(field) for field in nonblanks_filter]
        # we need to alter the query to look for blanks or limit to non-blanks
        # if specified
        
        # modify the query to include
        blanks_queries = []
        if blanks_filter:
            #blanks_queries = []
            for blank in blanks_filter:

                blanks_queries.append({"bool":
                                       {"should": [
                                           {"term": {blank + ".raw": ""}},
                                           {"missing": {"field": blank}}
                                       ]
                                       }
                                       })


        archived = request.GET.get("archived", None)
        if archived == "true":
            blanks_queries.append({"term": {"properties.archived": "true"}})
        else:
            blanks_queries.append({"bool":
                                   {"should": [{"term": {"properties.archived": "false"}},
                                               {"missing": {"field": "properties.archived"}}]}
                                   })


        nonblanks_queries = []
        if nonblanks_filter:

            for nonblank in nonblanks_filter:
                nonblanks_queries.append({"bool":
                                          {"should": [
                                              {"term": {
                                                  nonblank + ".raw": ""}},
                                              {"missing":
                                               {"field": nonblank}}
                                          ]
                                          }
                                          })
        modified_query = {
            "bool": {
                "must":  blanks_queries,
                "must_not": nonblanks_queries,
            },
        }
        if query:
            modified_query["bool"]["must"] += [query] 

        datestart = request.GET.get("created__gte", None)
        dateend = request.GET.get("created__lte", None)
        uoxs = request.GET.get("related_molregno__chembl__chembl_id__in", None)
        created_by = request.GET.get("created_by", "").split()
        if len(created_by) > 0:
            for id in created_by:
                modified_query["bool"]["must"] += [{
                                                        "term": {
                                                            "createdById" : int(id) 
                                                        }
                                                    }]

        projkeys = request.GET.get("project__project_key__in", "")
        if projkeys:
            crp = ChemregProjectResource()
            uri = crp.get_resource_uri()
            proj_ids = Project.objects.filter(project_key__in=projkeys.split(",")).values_list("pk", flat=True)
            pq = {"terms": {"project.raw": ["%s/%d" % (uri, pid) for pid in proj_ids]}}
            modified_query["bool"]["must"] += [pq,]

        if dateend or datestart:
            rq = {
                "range" : {
                    "created" : {
                        
                    }
                }
            }
            if dateend:
                dateend += " 23:59:59"
                rq["range"]["created"]["lte"] = dateend
            if datestart:
                rq["range"]["created"]["gte"] = datestart
            # rq["range"]["created"]["format"] = "yyyy/MM/dd HH:mm:ss"
            modified_query["bool"]["must"] += [rq,]

        if uoxs:
            #check for either chemblid field or blindedbatchid field
            tq = {"terms": {"chemblId.raw" : uoxs.split(",")}}
            tq2 = {"terms": {"blindedBatchId.raw" : uoxs.split(",")}}
            modified_query["bool"]["must"] += [{"bool": {"should" : [tq, tq2]}}]
        
        creator = request.GET.get("creator_uri", None)
        if creator:
            tq = {"terms": {"creator.raw": [cr for cr in creator.split(",")]}}
            modified_query["bool"]["must"] += [tq,]
        pids = self._meta.authorization.project_ids(request)
        pq = elasticsearch_client.get_project_uri_terms(pids)
        modified_query["bool"]["must"] += pq
        es_request = {
            "version": True,
            "from": get_data.get("offset", 0),
            "size": get_data.get("limit", 50),
            "query" : 
                { "bool": 
                    {"must":
                        [
                            {"filtered" : {"filter" : modified_query}}
                        ]
                    }
                },
            "sort": json.loads(get_data.get("sorts", '[{"id": {"order": "desc", "unmapped_type" : "long"}}]'))
        }

        prefix = request.GET.get("custom__field__startswith", -1)
        if prefix != -1:
            #Here the request is being used to get the custom field values
            custom_field = request.GET.get("custom_field", None)

            search_filter, aggs = elasticsearch_client.get_cf_aggregation(prefix, 'custom_field_list.aggregation', custom_field)
            modified_query["bool"]["must"] += [search_filter]
            es_request["aggs"] = aggs
            es_request["size"] = 0


        custom_fields__kv_any = request.GET.get("search_custom_fields__kv_any", "")
        if custom_fields__kv_any:
            custom_q = elasticsearch_client.get_custom_fields_query_from_string(custom_fields__kv_any)
            modified_query["bool"]["must"] += [custom_q]

        textsearch = request.GET.get("textsearch", None)
        if textsearch:
            #textsearch is a query not a filter
            es_request["query"]["bool"]["must"] += [{
                    "multi_match" :
                    { 
                        "type": "phrase_prefix", 
                        "fields": ["custom_field_list.value",] , 
                        "query" : textsearch 
                    }
                }]
        multiple_batch_id = request.GET.get("multiple_batch_id", None)
        if multiple_batch_id:
            modified_query["bool"]["must"] += [{"term": {"multiple_batch_id" : multiple_batch_id}}]

        #If only saved searches are to be included then filter for them, otherwise filter for all but them. Note that elasticsearch looks for 1 in a boolean data type field
        if kwargs.get("saved_search_projects_only", False):
            modified_query["bool"]["must"] += [{"term": {"projectfull.project_type.saved_search_project_type": "true" }}]
        else:
            modified_query["bool"]["must"] += [{"bool": {"must_not" :[{"term": {"projectfull.project_type.saved_search_project_type": "true" }}]}}]
        index = elasticsearch_client.get_main_index_name()
        es_serializer = CBHCompoundBatchElasticSearchSerializer()
        es_serializer.convert_query(es_request)
        bundledata = elasticsearch_client.get(index, es_request, {})
        if prefix != -1:
            data = [{"value": uox, "label": self.labelify_aggregate(
            uox, single_field=custom_field), } for uox in bundledata["aggregations"]]
            if custom_field:
                for item in data:
                    item["value"] = item["label"]
            serialized = json.dumps(data)
            rc = HttpResponse(
                content=serialized )
            return rc

        bundledata["objects"] = [
                    es_serializer.to_python_ready_data(d) for d in bundledata["objects"]
                ]
        bundledata = self.alter_list_data_to_serialize(request, bundledata)

        return self.create_response(request, bundledata)

    def get_cached_temporary_batch_data(self, multi_batch_id, get_data, request, bundledata={}):
        es_request = {
            "from": get_data.get("offset", 0),
            "size": get_data.get("limit", 50),
            "filter": json.loads(get_data.get("query", '{ "match_all" : {}}')),
            "sort": json.loads(get_data.get("sorts", '[{"id": {"order": "asc"}}]'))
        }
        index = elasticsearch_client.get_temp_index_name(
            request, multi_batch_id)
        es_serializer = CBHCompoundBatchElasticSearchSerializer()
        bundledata = elasticsearch_client.get(index, es_request, bundledata)
        bundledata["objects"] = [
            es_serializer.to_python_ready_data(d) for d in bundledata["objects"]]
        return bundledata

    def get_object_list(self, request):
        return super(CBHCompoundBatchResource, self).get_object_list(request)


def deepgetattr(obj, attr, ex):
    """Recurses through an attribute chain to get the ultimate value."""
    try:
        return reduce(getattr, attr.split('.'), obj)
    except:
        return ex

class CBHSavedSearchResource(CBHCompoundBatchResource):
    project = fields.ForeignKey(
        ChemregProjectResource, 'project', blank=False, null=False)

    class Meta(CBHCompoundBatchResource.Meta):
        resource_name = 'cbh_saved_search'
        es_index_name = "chemreg_chemical_index"


    def get_list_elasticsearch(self, request, **kwargs):
        return super(CBHSavedSearchResource, self).get_list_elasticsearch(request, saved_search_projects_only=True)



    


class CBHCompoundMultipleBatchResource(ModelResource):
    #comp_batch = fields.ForeignKey(CBHCompoundBatchResource, 'cbh_compound_batches', blank=False, null=False)
    #batches = fields.ToManyField(CBHCompoundBatchResource, 'batches', full=True)
    class Meta:
        filtering = {
            "created_by": ALL_WITH_RELATIONS,
            "project": ALL_WITH_RELATIONS,
        }
        always_return_data = True
        queryset = CBHCompoundMultipleBatch.objects.all()
        resource_name = 'cbh_multiple_batches'
        authorization = ProjectAuthorization()
        include_resource_uri = False
        allowed_methods = ['get']
        default_format = 'application/json'
        authentication = SessionAuthentication()

    def apply_filters(self, request, applicable_filters):
        pids = self._meta.authorization.project_ids(request)
        dataset = self.get_object_list(request).filter(
            **applicable_filters).filter(project_id__in=set(pids))

        return dataset.order_by("-created")

    def get_object_list(self, request):
        return super(CBHCompoundMultipleBatchResource, self).get_object_list(request).defer('uploaded_data').prefetch_related(Prefetch("project"))


class CBHCompoundBatchUpload(ModelResource):

    class Meta:
        excludes = ['uploaded_data']
        always_return_data = True
        queryset = FlowFile.objects.all()
        resource_name = 'cbh_batch_upload'
        authorization = Authorization()
        include_resource_uri = False
        allowed_methods = ['get', 'post', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/headers/?$" % self._meta.resource_name,
                self.wrap_view('return_headers'), name="api_compound_batch_headers"),
        ]

    def return_headers(self, request, **kwargs):
        deserialized = self.deserialize(request, request.body, format=request.META.get(
            'CONTENT_TYPE', 'application/json'))

        deserialized = self.alter_deserialized_detail_data(
            request, deserialized)
        bundle = self.build_bundle(
            data=dict_strip_unicode_keys(deserialized), request=request)
        request_json = bundle.data
        file_name = request_json['file_name']
        session_key = request.COOKIES[settings.SESSION_COOKIE_NAME]
        correct_file = self.get_object_list(request).get(
            identifier="%s-%s" % (session_key, file_name))

        header_json = {}
        if (correct_file.extension == ".sdf"):
            # read in the file
            headers = get_all_sdf_headers(correct_file.file.name)
        elif(correct_file.extension in (".xls", ".xlsx")):
            # read in excel file, use pandas to read the headers
            df = pd.read_excel(correct_file.file)
            headers = list(df)
        # this converts to json in preparation to be added to the response
        bundle.data["headers"] = list(set(headers))
        # send back
        # we should allow SD file uploads with no meta data
        if (len(headers) == 0 and correct_file.extension in (".xls", ".xlsx")):
            raise BadRequest("no_headers")
        return self.create_response(request, bundle, response_class=http.HttpAccepted)


def get_all_sdf_headers(filename):
    from subprocess import Popen, PIPE
    from shlex import split
    p1 = Popen(split('grep "^>" %s' % filename), stdout=PIPE)
    p2 = Popen(split('cut -d "<" -f2'), stdin=p1.stdout, stdout=PIPE)
    p3 = Popen(split('cut -d ">" -f1'), stdin=p2.stdout, stdout=PIPE)
    p4 = Popen(split('sort'), stdin=p3.stdout, stdout=PIPE)
    p5 = Popen(split('uniq'), stdin=p4.stdout, stdout=PIPE)
    out = p5.communicate()
    return [i for i in out[0].split("\n") if i]
