#!/usr/bin/python
# -*- coding: utf-8 -*-
from tastypie.resources import ModelResource, ALL_WITH_RELATIONS, ALL
from django.conf import settings
from django.conf.urls import url
from django.http import HttpResponse
from tastypie.http import  HttpConflict
from tastypie import fields
from cbh_core_model.models import Project, PinnedCustomField, CustomFieldConfig
from cbh_core_api.serializers import CustomFieldXLSSerializer
from cbh_core_api.resources import UserResource
from tastypie.authorization  import Authorization
from cbh_core_api.authorization import ProjectListAuthorization
from tastypie.authentication import SessionAuthentication
from tastypie.paginator import Paginator
from tastypie.exceptions import ImmediateHttpResponse
import json
import copy
import time
from django.core.urlresolvers import reverse
from cbh_core_api.serializers import CustomFieldsSerializer
from django.db.models import Prefetch
from cbh_core_api.resources import ProjectTypeResource, \
    CustomFieldConfigResource, UserHydrate
from django.contrib.auth.models import User
import six

def build_content_type(format, encoding='utf-8'):
    """
    Appends character encoding to the provided format if not already present.
    """

    if 'charset' in format:
        return format
    return '%s; charset=%s' % (format, encoding)


class ChemRegDataPointProjectFieldResource(ModelResource):

    """Provides the schema information about a field that is required by front end apps"""

    edit_form = fields.DictField(
        null=True, blank=False,  help_text=None)
    edit_schema = fields.DictField(
        null=True, blank=False,  help_text=None)
 
    class Meta:
        queryset = PinnedCustomField.objects.all()
        always_return_data = True
        resource_name = 'cbh_chemreg_datapoint_fields'
        include_resource_uri = True
        allowed_methods = ['get', 'post', 'patch', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        authorization = Authorization()
        level = None
        description = {'api_dispatch_detail' : '''
Provides information about the data types present in the flexible schema of the datapoint table
For each field a set of attributes are returned:

hide_form/schema - an angular schema form element that can be used to hide this column from view
edit_form /schema - an angular schema form element that can be used to edit this field 

assuming it is edited as part of a larger data form classification object
- To change the key of the json schema then change the get_namespace method

filter_form/schema - an angular schema form element that can be used to hide this filter this field

exclude_form /schema an angular schema form element that can be used to hide this exclude values from this field

sort_form /schema an angular schema form element that can be used to hide this exclude values from this field

Things still to be implemented:

actions form - would be used for mapping functions etc
autocomplete urls
        ''',

                       'api_dispatch_list' : '''
Provides information about the data types present in the flexible schema of the datapoint table
For each field a set of attributes are returned:

hide_form/schema - an angular schema form element that can be used to hide this column from view
edit_form /schema - an angular schema form element that can be used to edit this field 

assuming it is edited as part of a larger data form classification object
- To change the key of the json schema then change the get_namespace method

filter_form/schema - an angular schema form element that can be used to hide this filter this field

exclude_form /schema an angular schema form element that can be used to hide this exclude values from this field

sort_form /schema an angular schema form element that can be used to hide this exclude values from this field

Things still to be implemented:

actions form - would be used for mapping functions etc
autocomplete urls
        '''
         }



    def is_authenticated(self, request):
        """
        Handles checking if the user is authenticated and dealing with
        unauthenticated users.

        Mostly a hook, this uses class assigned to ``authentication`` from
        ``Resource._meta``.
        """
        # Authenticate the request as needed.
        return True


    def save(self, bundle, skip_errors=False):
        if bundle.via_uri:
            return bundle

        self.is_valid(bundle)

        if bundle.errors and not skip_errors:
            raise ImmediateHttpResponse(response=self.error_response(bundle.request, bundle.errors))

        # Check if they're authorized.
        # if bundle.obj.pk:
        #     self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        # else:
        #     self.authorized_create_detail(self.get_object_list(bundle.request), bundle)

        # Save FKs just in case.
        self.save_related(bundle)

        # Save the main object.
        obj_id = self.create_identifier(bundle.obj)
       
        if obj_id not in bundle.objects_saved or bundle.obj._state.adding:
            bundle.obj.save()
            bundle.objects_saved.add(obj_id)

        # Now pick up the M2M bits.
        m2m_bundle = self.hydrate_m2m(bundle)
        self.save_m2m(m2m_bundle)
        return bundle

    def get_schema(self, request, **kwargs):
        """
        Returns a serialized form of the schema of the resource.
        Calls ``build_schema`` to generate the data. This method only responds
        to HTTP GET.
        Should return a HttpResponse (200 OK).
        """
        # self.method_check(request, allowed=['get'])
        # self.is_authenticated(request)
        # self.throttle_check(request)
        # self.log_throttled_access(request)
        # bundle = self.build_bundle(request=request)
        # self.authorized_read_detail(self.get_object_list(bundle.request), bundle)
        return self.create_response(request, self.build_schema())

    def get_namespace(self, bundle):
        '''
            Hook to return the dotted path to this field based on the level and the name of the field
            The level name is formatted in the dehydrate method of the DataFormConfigResource
        '''
        return "{level}.project_data.%s" % (bundle.obj.get_space_replaced_name)

    def get_namespace_for_action_key(self, bundle, action_type):
        return action_type


    def dehydrate_edit_form(self, bundle):
        '''  Slightly different implementation of this        '''
        if bundle.request.GET.get("empty", False):
            return {}
        data = bundle.obj.field_values[1]
        data["key"] = bundle.obj.name
        if bundle.obj.UISELECTTAG in bundle.obj.field_type:
            data['options'] = {'refreshDelay': 0,
                                'async': {'url': "%s" % reverse('api_get_list_elasticsearch',
                                                         kwargs={'resource_name': 'cbh_compound_batches',
                                                                 'api_name': settings.WEBSERVICES_NAME})}
                                
                                }
        return {"form": [data]}

    def dehydrate_edit_schema(self, bundle):
        '''          '''
        if bundle.request.GET.get("empty", False):
            return {}
        return {"properties": {bundle.obj.name: bundle.obj.field_values[0]}}


    def authorized_update_detail(self, object_list, bundle):
        """
        Handles checking of permissions to see if the user has authorization
        to PUT this resource.
        """

        return True

    def authorized_create_detail(self, object_list, bundle):
        """
        Handles checking of permissions to see if the user has authorization
        to PUT this resource.
        """

        return True

class ChemRegCustomFieldConfigResource(UserHydrate, ModelResource):

    '''Return only the project type and custom field config name as returning the full field list would be '''
    data_type = fields.ForeignKey("cbh_core_api.resources.DataTypeResource",
                                  'data_type',  null=True, blank=False, default=None, full=True)
    project_data_fields = fields.ToManyField(ChemRegDataPointProjectFieldResource, 
        attribute="pinned_custom_field",null=True, blank=False, default=None, full=True)
    created_by = fields.ForeignKey(
        "cbh_core_api.resources.UserResource", 'created_by')

    class Meta:
        object_class = CustomFieldConfig
        queryset = CustomFieldConfig.objects.select_related(
            "created_by", "data_type",)
        excludes = ("schemaform")
        include_resource_uri = True
        resource_name = 'cbh_chemreg_custom_field_config'
        authentication = SessionAuthentication()
        authorization = Authorization()
        include_resource_uri = True
        default_format = 'application/json'
        serializer = CustomFieldXLSSerializer()
        # serializer = Serializer()
        filtering = {"id": ALL}
        allowed_methods = ['get', 'post', 'put', 'patch']
        description = {'api_dispatch_detail' : '''
Provides data about a single level of a data form config

data_type: A string to describe what "sort" of data this is (fields will generally be the same as other objects of this data type but that is up to the curator)
project_data_fields:
The fields that are in this particular custom field config:
    Provides information about the data types present in the flexible schema of the datapoint table
    For each field a set of attributes are returned:

    hide_form/schema - an angular schema form element that can be used to hide this column from view
    edit_form /schema - an angular schema form element that can be used to edit this field 

    assuming it is edited as part of a larger data form classification object
    - To change the key of the json schema then change the get_namespace method

    filter_form/schema - an angular schema form element that can be used to hide this filter this field
    
    exclude_form /schema an angular schema form element that can be used to hide this exclude values from this field
   
    sort_form /schema an angular schema form element that can be used to hide this exclude values from this field
   
   Things still to be implemented:

    actions form - would be used for mapping functions etc
    autocomplete urls
        ''',

                       'api_dispatch_list' : '''
Provides data about a single level of a data form config

data_type: A string to describe what "sort" of data this is (fields will generally be the same as other objects of this data type but that is up to the curator)
project_data_fields:
The fields that are in this particular custom field config:
    Provides information about the data types present in the flexible schema of the datapoint table
    For each field a set of attributes are returned:

    hide_form/schema - an angular schema form element that can be used to hide this column from view
    edit_form /schema - an angular schema form element that can be used to edit this field 

    assuming it is edited as part of a larger data form classification object
    - To change the key of the json schema then change the get_namespace method

    filter_form/schema - an angular schema form element that can be used to hide this filter this field
    
    exclude_form /schema an angular schema form element that can be used to hide this exclude values from this field
   
    sort_form /schema an angular schema form element that can be used to hide this exclude values from this field
   
   Things still to be implemented:

    actions form - would be used for mapping functions etc
    autocomplete urls
        '''
                       }


    def hydrate_project_data_fields(self, bundle):
        """Add the position number to the related bundle"""
        for index, item in enumerate(bundle.data["project_data_fields"]):
            if hasattr(item, 'obj'):
                item.obj.position = index
        return bundle


    def get_schema(self, request, **kwargs):
        """
        Returns a serialized form of the schema of the resource.
        Calls ``build_schema`` to generate the data. This method only responds
        to HTTP GET.
        Should return a HttpResponse (200 OK).
        """
        return self.create_response(request, self.build_schema())

    def create_response(self, request, data, response_class=HttpResponse, **response_kwargs):
        """
        Extracts the common "which-format/serialize/return-response" cycle.
        Mostly a useful shortcut/hook.
        """

        desired_format = self.determine_format(request)
        serialized = self.serialize(request, data, desired_format)
        rc = response_class(content=serialized, content_type=build_content_type(
            desired_format), **response_kwargs)

        if(desired_format == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'):
            rc['Content-Disposition'] = 'attachment; filename=project_data_explanation.xlsx'
        return rc





class ChemregProjectResource(UserHydrate, ModelResource):

    project_type = fields.ForeignKey(
        ProjectTypeResource, 'project_type', blank=False, null=False, full=True)
    custom_field_config = fields.ForeignKey(ChemRegCustomFieldConfigResource,
                                            'custom_field_config', blank=False, null=False, full=True)
    valid_cache_get_keys = ['format', 'limit', 'project_key',
                            'schemaform']
    assays_configured = fields.BooleanField(default=False)
    created_by = fields.ForeignKey(
        "cbh_core_api.resources.UserResource", 'created_by')

    class Meta:

        queryset = Project.objects.all()
        authentication = SessionAuthentication()
        paginator_class = Paginator
        allowed_methods = ['get', 'post', 'patch', 'put']
        resource_name = 'cbh_projects'
        authorization = ProjectListAuthorization()
        include_resource_uri = True
        default_format = 'application/json'
        serializer = CustomFieldsSerializer()
        filtering = {
            'project_key': ALL_WITH_RELATIONS,
            'project_type': ALL_WITH_RELATIONS,
        }
        always_return_data=True




    def dehydrate_assays_configured(self, bundle):
        return bundle.obj.enabled_forms.count() > 0

    def save_related(self, bundle):
        from django.db import IntegrityError
        try:
            return super(ChemregProjectResource, self).save_related(bundle)
        except IntegrityError, e:
            raise ImmediateHttpResponse(HttpConflict("Project with that name already exists"))




    def post_list(self, request, **kwargs):
        from django.db import IntegrityError
        try:
            return super(ChemregProjectResource, self).post_list(request, **kwargs)
        except IntegrityError, e:
            raise ImmediateHttpResponse(HttpConflict("Project with that name already exists"))


    def patch_detail(self, request, **kwargs):
        from django.db import IntegrityError
        try:
            return super(ChemregProjectResource, self).patch_detail(request, **kwargs)
        except IntegrityError, e:
            raise ImmediateHttpResponse(HttpConflict("Project with that name already exists"))



    def patch_list(self, request, **kwargs):
        from django.db import IntegrityError
        try:
            return super(ChemregProjectResource, self).patch_list(request, **kwargs)
        except IntegrityError, e:
            raise ImmediateHttpResponse(HttpConflict("Project with that name already exists"))





    def get_object_list(self, request):
        return super(ChemregProjectResource,
                     self).get_object_list(request).prefetch_related(Prefetch('project_type'
                                                                              )).order_by('-modified')

    def prepend_urls(self):
        return [url(r"^(?P<resource_name>%s)/custom_fields/?$"
                    % self._meta.resource_name,
                    self.wrap_view('get_custom_fields'),
                    name='get_custom_fields')]

    def get_custom_fields(self, request):
        return super(ChemregProjectResource,
                     self).get_object_list(request).prefetch_related(Prefetch('custom_field_config'
                                                                              ))

    def get_searchform(self, bundle):
        '''Note that the form here was expected to have the UOx id as the first item
           PB edit - we are now looking it up by key in front end - UOX position no longer needs to be fixed'''
        ur = UserResource()
        uri = ur.get_resource_uri()
        return {
            'cf_form': [{
                'htmlClass': 'col-sm-10',
                'key': 'search_custom_fields__kv_any',
                'disableSuccessState': True,
                'feedback': False,
                'options': {'refreshDelay': 0,
                            'async': {'url': reverse('api_get_elasticsearch_autocomplete',
                                                     kwargs={'resource_name': 'cbh_compound_batches',
                                                             'api_name': settings.WEBSERVICES_NAME})}},
            }],
            'cf_schema': {'required': [], 'type': 'object',
                          'properties': {'search_custom_fields__kv_any': {
                              'type': 'array',
                              'format': 'uiselect',
                              'items': [],
                              'placeholder': 'Filter project data',
                              'title': 'Project data values:',
                          }}},
            #add a simple form for the save search interface to allow uiselect of projects and groups
            'savesearch_form': [
                # {
                #     'key': 'project__project_key__in',
                #     'placeholder': 'Select project',
                #     'htmlClass': 'col-xs-12',
                #     'feedback': False,
                #     #'description': 'Search for projects in order to limit the choice of fields on show. Select a single project if you want to edit data.',
                #     'disableSuccessState': True,
                #     #'validationMessage': {'default': 'Please select a project if you wish to edit data.'}
                # },
                {
                    'key': 'alias',
                    'placeholder': 'e.g. My Saved Search',
                    'htmlClass': 'col-xs-12',
                    'validationMessage': 'Please add an alias for this search',
                },
            ],
            'savesearch_schema': {
                          
                          'type': 'object',
                          'properties': {
                                'alias': {
                                    'title': 'Save search as...',
                                    'type': 'string',
                                },
                            },
                            'required': ['alias'], 
                        },
            'simple_form': [

                {
                    'key': 'creator_uri',
                    'htmlClass': 'col-md-4 col-xs-6',
                    'placeholder': 'Select users to search',
                    'feedback': False,


                },
                {
                    'key': 'project__project_key__in',
                    'placeholder': 'Select projects to search',
                    'htmlClass': 'col-md-4 col-xs-6',
                    'feedback': False,
                    'description': 'Search for projects in order to limit the choice of fields on show. Select a single project if you want to edit data.',
                    'disableSuccessState': True,
                    'validationMessage': {'default': 'Please select a project if you wish to edit data.'}
                },
                {
                    'htmlClass': 'col-md-4 col-xs-6',
                    'key': 'search_custom_fields__kv_any',
                    'disableSuccessState': True,
                    'help': 'Searching using this filter will bring back results that match an OR pattern within the same data category, with AND across data categories, i.e. results which contain this item within category a OR that item within category a AND that item within category b.',
                    'feedback': False,
                    'options': {'refreshDelay': 0,
                                'async': {'url': reverse('api_get_list_elasticsearch',
                                                         kwargs={'resource_name': 'cbh_compound_batches',
                                                                 'api_name': settings.WEBSERVICES_NAME})},
                                
                                },
                },
            ],
            'simple_schema': {
                          'required': [], 'type': 'object',
                          'properties': {
                                'search_custom_fields__kv_any': {
                                  'type': 'array',
                                  'format': 'uiselect',
                                  'items': [],
                                  'placeholder': 'Filter project data',
                                  'title': 'Project data values:',
                                },
                                'creator_uri': {
                                    'type': 'array',
                                    'format': 'uiselect',
                                       'title': 'Compound batch created by',
                                        'type': 'array',
                                        'format': 'uiselect',
                                        'htmlClass': 'col-md-6 col-xs-6',
                                        'placeholder': 'Search user who created the batch',
                                        'options': {'searchDescriptions': False},
                                        'items':  sorted([
                                            {'label': user.first_name + " " + user.last_name , "value" : uri + '/' + str(user.id) } if user.first_name  else {'label': user.username , "value" : uri + '/' + str(user.id) }
                                            for user in User.objects.exclude(pk=-1)
                                        ], key=lambda k: k['label'].lower())
                                },
                                'project__project_key__in': {
                                    'title': 'Project',
                                    'type': 'array',
                                    'format': 'uiselect',
                                    'items': [{'label': p.obj.name,
                                               'value': p.obj.project_key} for p in
                                              bundle['objects']],
                                },
                          }},
            
            'form': [
                {
                    'key': 'related_molregno__chembl__chembl_id__in',
                    'title': '%s ID' % settings.ID_PREFIX,
                    'placeholder': 'Search multiple IDs',
                    'feedback': False,
                        'htmlClass': 'col-md-6 col-xs-6',

                    'options': {'refreshDelay': 0,
                                'async': {'url': reverse('api_get_elasticsearch_ids',
                                                         kwargs={'resource_name': 'cbh_compound_batches',
                                                                 'api_name': settings.WEBSERVICES_NAME})}},
                },
                {
                    'key': 'multiple_batch_id',
                    'htmlClass': 'col-md-6 col-xs-6',
                    'disableSuccessState': True,
                    'feedback': False,
                },
                {
                    'key': 'dateStart',
                    'type': 'datepicker',
                    'minDate': '2004-01-01',
                    'htmlClass': 'col-md-6 col-xs-6',
                    'disableSuccessState': True,
                    'feedback': False,
                    'pickadate': {
                        'selectYears': True,
                        'selectMonths': True
                        },
                },
                {
                    'key': 'dateEnd',
                    'type': 'datepicker',
                    'minDate': '2004-01-01',
                        'htmlClass': 'col-md-6 col-xs-6',
                    'disableSuccessState': True,
                    'feedback': False,
                    'pickadate': {'selectYears': True,
                                  'selectMonths': True},
                }, 
                {
                    'htmlClass': 'col-md-6 col-xs-6',
                    'disableSuccessState': True,
                    'feedback': False,
                    'key': 'functional_group',

                },
                {
                    'key': 'smiles',
                    'placeholder': 'Search SMILES or SMARTS string',
                    'append': 'today',
                    'feedback': False,
                    'htmlClass': 'col-md-6 col-xs-6',
                    'disableSuccessState': True,
                },
                {
                    'key': 'substruc',
                    'style': {'selected': 'btn-success',
                              'unselected': 'btn-default'},
                        'htmlClass': 'col-md-6 col-xs-6',
                    'type': 'radiobuttons',
                    'disableSuccessState': True,
                    'feedback': False,
                    'titleMap': [{'value': 'with_substructure',
                                  'name': 'Substructure'},
                                 {'value': 'flexmatch',
                                  'name': 'Exact Match'}],
                },
                {
                    'key': 'archived',
                    'style': {'selected': 'btn-success',
                              'unselected': 'btn-default'},
                        'htmlClass': 'col-md-6 col-xs-6',
                    'type': 'radiobuttons',
                    'disableSuccessState': True,
                    'feedback': False,
                    'titleMap': [
                        {'value': 'false',
                                  'name': 'Normal mode'},
                        {'value': 'true',
                                  'name': 'Archive mode'},]
                                 
                },
            ],
            'schema': {'required': [], 'type': 'object', 'properties': {
                'related_molregno__chembl__chembl_id__in': {
                    'type': 'array',
                    'format': 'uiselect',
                    
                },

          
                
                'multiple_batch_id': {'title': 'Upload ID',
                                      'type': 'string'},

                'functional_group': {
                    'title': 'Functional Group',
                    'type': 'string',
                    'format': 'uiselect',
                    'placeholder': 'Search chemical groups',
                    'options': {'searchDescriptions': False},
                    'default': '',
                    'copyValueTo': 'smiles',
                    'items': [{'label': 'None', 'value': ''}] + sorted([
                        {'label': 'Alkyl Carbon', 'value': '[CX4]'},
                        {'label': 'Allenic Carbon',
                         'value': '[$([CX2](=C)=C)]'},
                        {'label': 'Vinylic Carbon',
                         'value': '[$([CX3]=[CX3])]'},
                        {'label': 'Acetylenic Carbon',
                         'value': '[$([CX2]#C)]'},
                        {'label': 'Arene', 'value': 'c'},
                        {'label': 'Carbonyl group. Low specificity',
                         'value': '[CX3]=[OX1]'},
                        {'label': 'Carbonyl group',
                         'value': '[$([CX3]=[OX1]),$([CX3+]-[OX1-])]'},
                        {'label': 'Carbonyl with Carbon',
                         'value': '[CX3](=[OX1])C'},
                        {'label': 'Carbonyl with Nitrogen.',
                         'value': '[OX1]=CN'},
                        {'label': 'Carbonyl with Oxygen.',
                         'value': '[CX3](=[OX1])O'},
                        {'label': 'Acyl Halide',
                         'value': '[CX3](=[OX1])[F,Cl,Br,I]'},
                        {'label': 'Aldehyde', 'value': '[CX3H1](=O)[#6]'
                         },
                        {'label': 'Anhydride',
                         'value': '[CX3](=[OX1])[OX2][CX3](=[OX1])'},
                        {'label': 'Amide',
                         'value': '[NX3][CX3](=[OX1])[#6]'},
                        {'label': 'Amidinium',
                         'value': '[NX3][CX3]=[NX3+]'},
                        {'label': 'Carbamate.',
                         'value': '[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]'},
                        {'label': 'Carbamic ester',
                         'value': '[NX3][CX3](=[OX1])[OX2H0]'},
                        {'label': 'Carbamic acid.',
                         'value': '[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]'
                         },
                        {'label': 'Carboxylate Ion.',
                         'value': '[CX3](=O)[O-]'},
                        {'label': 'Carbonic Acid or Carbonic Ester',
                         'value': '[CX3](=[OX1])(O)O'},
                        {'label': 'Carbonic Acid or Carbonic Acid-Ester', 'value': '[CX3](=[OX1])([OX2])[OX2H,OX1H0-1]'
                         },
                        {'label': 'Carbonic Ester (carbonic acid diester)',
                         'value': 'C[OX2][CX3](=[OX1])[OX2]C'},
                        {'label': 'Carboxylic acid',
                         'value': '[CX3](=O)[OX2H1]'},
                        {'label': 'Carboxylic acid or conjugate base.',
                         'value': '[CX3](=O)[OX1H0-,OX2H1]'},
                        {'label': 'Cyanamide',
                         'value': '[NX3][CX2]#[NX1]'},
                        {'label': 'Ester Also hits anhydrides',
                         'value': '[#6][CX3](=O)[OX2H0][#6]'},
                        {'label': 'Ketone', 'value': '[#6][CX3](=O)[#6]'
                         },
                        {'label': 'Ether', 'value': '[OD2]([#6])[#6]'},
                        {'label': 'Hydrogen Atom', 'value': '[H]'},
                        {'label': 'Not a Hydrogen Atom',
                         'value': '[!#1]'},
                        {'label': 'Proton', 'value': '[H+]'},
                        {'label': 'Mono-Hydrogenated Cation',
                         'value': '[+H]'},
                        {'label': 'Not Mono-Hydrogenated',
                         'value': '[!H] or [!H1]'},
                        {'label': 'Primary or secondary amine, not amide.',
                            'value': '[NX3;H2,H1;!$(NC=O)]'},
                        {'label': 'Enamine', 'value': '[NX3][CX3]=[CX3]'
                         },
                        {'label': 'Primary amine, not amide.',
                         'value': "[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6] Not amide (C not double bonded to a hetero-atom), not ammonium ion (N must be 3-connected), not ammonia (N's H-count can't be 3), not cyanamide (C not triple bonded to a hetero-atom)"
                         },
                        {'label': 'Two primary or secondary amines',
                         'value': '[NX3;H2,H1;!$(NC=O)].[NX3;H2,H1;!$(NC=O)]'
                         },
                        {'label': 'Enamine or Aniline Nitrogen',
                         'value': '[NX3][$(C=C),$(cc)]'},
                        {'label': 'Azide group.',
                         'value': '[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]'
                         },
                        {'label': 'Azide ion.',
                         'value': '[$([NX1-]=[NX2+]=[NX1-]),$([NX1]#[NX2+]-[NX1-2])]'
                         },
                        {'label': 'Nitrogen.', 'value': '[#7]'},
                        {'label': 'Azo Nitrogen. Low specificity.',
                         'value': '[NX2]=N'},
                        {'label': 'Azo Nitrogen.diazene',
                         'value': '[NX2]=[NX2]'},
                        {'label': 'Azoxy Nitrogen.',
                         'value': '[$([NX2]=[NX3+]([O-])[#6]),$([NX2]=[NX3+0](=[O])[#6])]'
                         },
                        {'label': 'Diazo Nitrogen',
                         'value': '[$([#6]=[N+]=[N-]),$([#6-]-[N+]#[N])]'
                         },
                        {'label': 'Azole.',
                         'value': '[$([nr5]:[nr5,or5,sr5]),$([nr5]:[cr5]:[nr5,or5,sr5])]'
                         },
                        {'label': 'Hydrazine H2NNH2',
                         'value': '[NX3][NX3]'},
                        {'label': 'Hydrazone C=NNH2',
                         'value': '[NX3][NX2]=[*]'},
                        {'label': 'Substituted imine',
                         'value': '[CX3;$([C]([#6])[#6]),$([CH][#6])]=[NX2][#6]'
                         },
                        {'label': 'Substituted or un-substituted imine',
                         'value': '[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]'
                         },
                        {'label': 'Iminium', 'value': '[NX3+]=[CX3]'},
                        {'label': 'Unsubstituted dicarboximide',
                         'value': '[CX3](=[OX1])[NX3H][CX3](=[OX1])'},
                        {'label': 'Substituted dicarboximide',
                         'value': '[CX3](=[OX1])[NX3H0]([#6])[CX3](=[OX1])'
                         },
                        {'label': 'Dicarboxdiimide',
                         'value': '[CX3](=[OX1])[NX3H0]([NX3H0]([CX3](=[OX1]))[CX3](=[OX1]))[CX3](=[OX1])'
                         },
                        {'label': 'Nitrate group',
                         'value': '[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]'
                         },
                        {'label': 'Nitrate Anion',
                         'value': '[$([OX1]=[NX3](=[OX1])[OX1-]),$([OX1]=[NX3+]([OX1-])[OX1-])]'
                         },
                        {'label': 'Nitrile', 'value': '[NX1]#[CX2]'},
                        {'label': 'Isonitrile', 'value': '[CX1-]#[NX2+]'
                         },
                        {'label': 'Nitro group.',
                         'value': '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8] Hits both forms.'
                         },
                        {'label': 'Two Nitro groups',
                         'value': '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8].[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]'
                         },
                        {'label': 'Nitroso-group',
                         'value': '[NX2]=[OX1]'},
                        {'label': 'N-Oxide',
                         'value': '[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]'
                         },
                        {'label': 'Hydroxyl', 'value': '[OX2H]'},
                        {'label': 'Hydroxyl in Alcohol',
                         'value': '[#6][OX2H]'},
                        {'label': 'Hydroxyl in Carboxylic Acid',
                         'value': '[OX2H][CX3]=[OX1]'},
                        {'label': 'Hydroxyl in H-O-P-',
                         'value': '[OX2H]P'},
                        {'label': 'Enol', 'value': '[OX2H][#6X3]=[#6]'
                         },
                        {'label': 'Phenol', 'value': '[OX2H][cX3]:[c]'
                         },
                        {'label': 'Enol or Phenol',
                         'value': '[OX2H][$(C=C),$(cc)]'},
                        {'label': 'Hydroxyl_acidic',
                         'value': '[$([OH]-*=[!#6])]'},
                        {'label': 'Peroxide groups.',
                         'value': '[OX2,OX1-][OX2,OX1-]'},
                        {'label': 'Phosphoric_acid groups.',
                         'value': '[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]'
                         },
                        {'label': 'Phosphoric_ester groups.',
                         'value': '[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]'
                         },
                        {'label': 'Carbo-Thiocarboxylate',
                         'value': '[S-][CX3](=S)[#6]'},
                        {'label': 'Carbo-Thioester',
                         'value': 'S([#6])[CX3](=O)[#6]'},
                        {'label': 'Thio analog of carbonyl',
                         'value': '[#6X3](=[SX1])([!N])[!N]'},
                        {'label': 'Thiol, Sulfide or Disulfide Sulfur',
                         'value': '[SX2]'},
                        {'label': 'Thiol', 'value': '[#16X2H]'},
                        {'label': 'Sulfur with at-least one hydrogen.',
                         'value': '[#16!H0]'},
                        {'label': 'Thioamide',
                         'value': '[NX3][CX3]=[SX1]'},
                        {'label': 'Sulfide', 'value': '[#16X2H0]'},
                        {'label': 'Mono-sulfide',
                         'value': '[#16X2H0][!#16]'},
                        {'label': 'Di-sulfide',
                         'value': '[#16X2H0][#16X2H0]'},
                        {'label': 'Two Sulfides',
                         'value': '[#16X2H0][!#16].[#16X2H0][!#16]'},
                        {'label': 'Sulfinate',
                         'value': '[$([#16X3](=[OX1])[OX2H0]),$([#16X3+]([OX1-])[OX2H0])]'
                         },
                        {'label': 'Sulfinic Acid',
                         'value': '[$([#16X3](=[OX1])[OX2H,OX1H0-]),$([#16X3+]([OX1-])[OX2H,OX1H0-])]'
                         },
                        {'label': 'Sulfone. Low specificity.',
                         'value': '[$([#16X4](=[OX1])=[OX1]),$([#16X4+2]([OX1-])[OX1-])]'
                         },
                        {'label': 'Sulfone. High specificity.',
                         'value': '[$([#16X4](=[OX1])(=[OX1])([#6])[#6]),$([#16X4+2]([OX1-])([OX1-])([#6])[#6])]'
                         },
                        {'label': 'Sulfonic acid. High specificity.',
                         'value': '[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]'
                         },
                        {'label': 'Sulfonate',
                         'value': '[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]'
                         },
                        {'label': 'Sulfonamide.',
                         'value': '[$([#16X4]([NX3])(=[OX1])(=[OX1])[#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[#6])]'
                         },
                        {'label': 'Carbo-azosulfone',
                         'value': '[SX4](C)(C)(=O)=N'},
                        {'label': 'Sulfonamide',
                         'value': '[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]'
                         },
                        {'label': 'Sulfoxide Low specificity.',
                         'value': '[$([#16X3]=[OX1]),$([#16X3+][OX1-])]'
                         },
                        {'label': 'Sulfoxide High specificity',
                         'value': '[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]'
                         },
                        {'label': 'Sulfate',
                         'value': '[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]'
                         },
                        {'label': 'Sulfuric acid ester (sulfate ester) Low specificity.',
                         'value': '[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]'
                         },
                        {'label': 'Sulfuric Acid Diester.',
                         'value': '[$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6]),$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6])]'
                         },
                        {'label': 'Sulfamate.',
                         'value': '[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2][#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2][#6])]'
                         },
                        {'label': 'Sulfamic Acid.',
                         'value': '[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2H,OX1H0-]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2H,OX1H0-])]'
                         },
                        {'label': 'Sulfenic acid.',
                         'value': '[#16X2][OX2H,OX1H0-]'},
                        {'label': 'Sulfenate.',
                         'value': '[#16X2][OX2H0]'},
                        {'label': 'Any carbon attached to any halogen',
                         'value': '[#6][F,Cl,Br,I]'},
                        {'label': 'Halogen', 'value': '[F,Cl,Br,I]'},
                        {'label': 'Sulfide', 'value': '[#16X2H0]'},
                        {'label': 'Mono-sulfide',
                         'value': '[#16X2H0][!#16]'},
                        {'label': 'Di-sulfide',
                         'value': '[#16X2H0][#16X2H0]'},
                        {'label': 'Hydrogen-bond acceptor',
                         'value': '[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]'
                         },
                        {'label': 'Hydrogen-bond donor.',
                         'value': '[!$([#6,H0,-,-2,-3])]'},
                    ], key=lambda k: k['label']),
                },
                'dateStart': {
                    'title': 'Added after',
                    'type': 'string',
                    'format': 'date',
                    'style': {'margin-right': '30px;'},
                },
                'dateEnd': {'title': 'Added before', 'type': 'string',
                            'format': 'date'},
                'smiles': {'title': 'SMILES or SMARTS',
                           'type': 'string'},
                'substruc': {
                    'title': 'Structural search type',
                    'type': 'string',
                    'enum': ['with_substructure', 'flexmatch'],
                    'default': 'with_substructure',
                },
                'search_custom_fields__kv_any': {
                    'type': 'array',
                    'format': 'uiselect',
                    'items': [],
                    'placeholder': 'Choose column and value...',
                    'title': 'Filter by project data values:',
                },
                'archived': {
                    'title': 'Enable/disable archive mode',
                    'type': 'string',
                    'enum': ['true', 'false'],
                    'default': 'false',
                },

            }},
        }

    def alter_list_data_to_serialize(self, request, bundle):
        '''Here we append a list of tags to the data of the GET request if the
        search fields are required'''

        userres = UserResource()
        userbundle = userres.build_bundle(obj=request.user,
                                          request=request)
        userbundle = userres.full_dehydrate(userbundle)
        bundle['user'] = userbundle.data
        self._meta.authorization.alter_project_data_for_permissions(bundle, request)

        if request.GET.get('schemaform', None):
            searchfields = set([])
            searchfield_items = []
            bundle['searchform'] = self.get_searchform(bundle, )
        return bundle


    def alter_detail_data_to_serialize(self, request, bundle):
        userres = UserResource()
        userbundle = userres.build_bundle(obj=request.user,
                                          request=request)
        userbundle = userres.full_dehydrate(userbundle)
        bundle.data['user'] = userbundle.data

        self._meta.authorization.alter_project_data_for_permissions(bundle, request)
                          
        return bundle

   



    def create_response(
        self,
        request,
        data,
        response_class=HttpResponse,
        **response_kwargs ):
        """
        Extracts the common "which-format/serialize/return-response" cycle.
        Mostly a useful shortcut/hook.
        """

        desired_format = self.determine_format(request)
        serialized = self.serialize(request, data, desired_format)
        rc = response_class(content=serialized,
                            content_type=build_content_type(desired_format),
                            **response_kwargs)

        if desired_format \
                == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet':
            rc['Content-Disposition'] = \
                'attachment; filename=project_data_explanation.xlsx'
        return rc



class NoCustomFieldsChemregProjectResource(ChemregProjectResource):
    custom_field_config = fields.ForeignKey(ChemRegCustomFieldConfigResource,
                                            'custom_field_config', blank=False, null=False)
    class Meta(ChemregProjectResource.Meta):
        pass