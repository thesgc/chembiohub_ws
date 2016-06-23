"""
Main interface for requesting the project JSON data from ChemBio Hub Platform
Also schema and configuration APIs
"""
#FlowFile relocation
import os
import datetime
from cbh_core_api.flowjs_settings import FLOWJS_PATH, FLOWJS_EXPIRATION_DAYS


import logging

# Get an instance of a logger
logger = logging.getLogger(__name__)

from django.middleware.csrf import get_token
from django.test import RequestFactory
from tastypie.resources import ALL_WITH_RELATIONS, ALL
from tastypie.resources import ModelResource
from tastypie.exceptions import BadRequest
from django.conf import settings
from django.conf.urls import *
from django.http import HttpResponse, QueryDict
import shortuuid
from tastypie.resources import ModelResource
from tastypie import fields
from django.contrib.auth.tokens import default_token_generator
from tastypie.exceptions import ImmediateHttpResponse

from django.contrib.auth.forms import PasswordResetForm
from cbh_core_model.models import CustomFieldConfig
from cbh_core_model.models import Project
from cbh_core_model.models import ProjectType
from cbh_core_model.models import SkinningConfig
from cbh_core_model.models import Invitation
from cbh_core_model.models import PinnedCustomField


from cbh_core_api.authorization import ProjectListAuthorization, InviteAuthorization, viewer_projects, ProjectPermissionAuthorization
from tastypie.authentication import SessionAuthentication
from tastypie.paginator import Paginator

from django.db.models import Prefetch

from tastypie.resources import ALL_WITH_RELATIONS
from tastypie.utils.mime import build_content_type
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.conf import settings
from cbh_core_api.authorization import get_all_project_ids_for_user
from django.conf import settings
from django.views.generic import FormView, View
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth import login as auth_login, logout as auth_logout
from tastypie.resources import ModelResource
from tastypie.authorization import Authorization
from tastypie import http
from django.contrib.auth.views import password_reset
from django.db import IntegrityError
try:
    # django >= 1.7
    from django.apps import apps
    get_model = apps.get_model
except ImportError:
    # django < 1.7
    from django.db.models import get_model

# If ``csrf_exempt`` isn't present, stub it.
try:
    from django.views.decorators.csrf import csrf_exempt
except ImportError:
    def csrf_exempt(func):
        return func
from  BeautifulSoup import BeautifulSoup
try:
    import defusedxml.lxml as lxml
except ImportError:
    lxml = None

try:
    WS_DEBUG = settings.WS_DEBUG
except AttributeError:
    WS_DEBUG = False


from tastypie.authentication import SessionAuthentication

from django.contrib.auth import get_user_model
import inflection
import six
import importlib

from django.views.generic import TemplateView
from django.views.decorators.csrf import csrf_exempt
from django.utils.decorators import method_decorator
from django.contrib.auth.forms import PasswordResetForm, loader, get_current_site, urlsafe_base64_encode, force_bytes

from urllib import urlencode
from django.core.mail import EmailMessage
from django.contrib.auth.models import Permission
import re
from django_hstore import hstore
#!/usr/bin/python
# -*- coding: utf-8 -*-
from tastypie.resources import ModelResource, ALL_WITH_RELATIONS, ALL
from django.conf import settings
from django.conf.urls import url
from django.middleware import csrf
from django.http import HttpResponse
from tastypie.http import  HttpConflict
from tastypie import fields
from cbh_core_model.models import Project, PinnedCustomField, CustomFieldConfig
from cbh_core_api.serializers import CustomFieldXLSSerializer
from tastypie.authorization  import Authorization
from cbh_core_api.authorization import ProjectListAuthorization, ProjectAuthorization
from tastypie.authentication import SessionAuthentication
from tastypie.paginator import Paginator
from tastypie.exceptions import ImmediateHttpResponse
import json
import copy
import time
from django.core.urlresolvers import reverse
from django.db.models import Prefetch
from django.contrib.auth.models import User
import six
from tastypie.serializers import Serializer
from cbh_core_api.tasks import remove_session_cached_projectlists
from django.contrib.auth.models import User

from cbh_utils import elasticsearch_client

def build_content_type(format, encoding='utf-8'):
    """
    Appends character encoding to the provided format if not already present.
    """

    if 'charset' in format:
        return format
    return '%s; charset=%s' % (format, encoding)


class UserHydrate(object):
    """Abstract class to make sure each 
    resource object can hydrate the user object against it"""
    def hydrate_created_by(self, bundle):
        User = get_user_model()

        if bundle.obj.id:
            pass
        else:
            user = User.objects.get(pk=bundle.request.user.pk)
            bundle.obj.created_by = user
        return bundle



class ChemRegDataPointProjectFieldResource(ModelResource):

    """Provides the schema information about a field that is required by front end apps"""

    edit_form = fields.DictField(
        null=True, blank=False,  help_text="an angular schema form element that can be used to edit this field ")
    edit_schema = fields.DictField(
        null=True, blank=False,   help_text="an angular schema form schema that can be used to edit this field ")
    view_form = fields.DictField(
        null=True, blank=False, help_text="possibly deprecated or never used"
        )
    handsontable = fields.DictField(
        null=True, blank=False, help_text="A JSON object which is combined with the other fields on the system to make a handsontable columns object"
        )
    

 
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

edit_form /schema - an angular schema form element that can be used to edit this field 



        ''',

                       'api_dispatch_list' : '''
Provides information about the data types present in the flexible schema of the datapoint table
For each field a set of attributes are returned:

edit_form /schema - an angular schema form element that can be used to edit this field 


        '''
         }

    def dehydrate_handsontable(self, bundle):
        """Generate the JSON object required by the handontable API"""
        return {
            "data" : "custom_fields.%s" % bundle.obj.name,
            "knownBy" : bundle.obj.name,
            "export_name" : bundle.obj.name,
            "renderer_named" : "customFieldRenderer",
            "className" : "htCenter htMiddle ",
            "editable" : True, # flags the field as editable so edit links appear in it
            "editor" : False    # ensures that the handsontable in cell editor is disabled

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
        """
        Cannot remember why the save function is "re implemented" here
        """
        if bundle.via_uri:
            return bundle

        self.is_valid(bundle)

        if bundle.errors and not skip_errors:
            raise ImmediateHttpResponse(response=self.error_response(bundle.request, bundle.errors))


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

        return self.create_response(request, self.build_schema())



    def dehydrate_edit_form(self, bundle):
        ''' Return the part of the edit form related to this field -
        The edit forms will be combined on the front end see 
        https://github.com/thesgc/ng-chem/blob/master/app/scripts/config.js'''
        if bundle.request.GET.get("empty", False):
            return {}
        data = bundle.obj.field_values[1]
        
        
        return {"form": [data]}

    def dehydrate_edit_schema(self, bundle):
        ''' Return the part of the edit schema related to this field -
        The edit schemas will be combined on the front end see 
        https://github.com/thesgc/ng-chem/blob/master/app/scripts/config.js'''
        if bundle.request.GET.get("empty", False):
            return {}
        return {"properties": {bundle.obj.name: bundle.obj.field_values[0]}}

    def dehydrate_view_form(self, bundle):
        '''possibly deprecated  '''
        if bundle.request.GET.get("empty", False):
            return {}
        data = bundle.obj.field_values[2]
        return {"form" : data}





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
    project_data_fields = fields.ToManyField(ChemRegDataPointProjectFieldResource, 
        attribute="pinned_custom_field",null=True, blank=False, default=None, full=True, help_text="List of the fields related to this cusotm field config, to be combined together to produce an edit schema")
    created_by = fields.ForeignKey(
        "cbh_core_api.resources.UserResource", 'created_by', help_text="The user who added this object")

    class Meta:
        object_class = CustomFieldConfig
        queryset = CustomFieldConfig.objects.select_related(
            "created_by",)
        excludes = ("schemaform")
        include_resource_uri = True
        resource_name = 'cbh_chemreg_custom_field_config'
        authentication = SessionAuthentication()
        authorization = Authorization()
        include_resource_uri = True
        default_format = 'application/json'
        serializer = CustomFieldXLSSerializer()
        filtering = {"id": ALL}
        allowed_methods = ['get', 'post', 'put', 'patch']
        description = {'api_dispatch_detail' : '''

        ''',

                       'api_dispatch_list' : '''

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
        "cbh_core_api.resources.ProjectTypeResource", 'project_type',  blank=False, null=False, full=True, help_text="A tag for the type of data this project stores")
    custom_field_config = fields.ForeignKey(ChemRegCustomFieldConfigResource,
                                            'custom_field_config', blank=False, null=False, full=True, help_text="The single custom field config object attached to this project")
    valid_cache_get_keys = ['format', 'limit', 'project_key',
                            'schemaform']
    created_by = fields.ForeignKey(
        "cbh_core_api.resources.UserResource", 'created_by', help_text="The user who created this Project")
    users_restricted_fields = fields.ListField(default=[], help_text="Possibly deprecated, was meant to list the restricted fields")
    flowjs_upload_url = fields.CharField(default="", help_text="The URL to be used when uploading data or attachments associated with this project")

    class Meta:

        queryset = Project.objects.all()
        authentication = SessionAuthentication()
        paginator_class = Paginator
        allowed_methods = ['get', 'post', 'patch', 'put']
        resource_name = 'cbh_projects'
        authorization = ProjectListAuthorization()
        include_resource_uri = True
        default_format = 'application/json'
        serializer = Serializer()
        filtering = {
            'project_key': ALL_WITH_RELATIONS,
            'project_type': ALL_WITH_RELATIONS,
            'id': ALL_WITH_RELATIONS,
        }
        always_return_data=True

    def dehydrate_flowjs_upload_url(self, bundle):
        try:
            return reverse('flowv2_upload', kwargs={'project_id': bundle.obj.id })
        except Exception:
            return ""



    def save_related(self, bundle):
        """Ensure that when saving the custom field config, any integrity errors from repeated names are passed up"""
        from django.db import IntegrityError
        try:
            return super(ChemregProjectResource, self).save_related(bundle)
        except IntegrityError, e:
            raise ImmediateHttpResponse(HttpConflict("Project with that name already exists"))


    def get_list(self, request, **kwargs):
        """
        Request the list response and if asked to, use the cache
        Cache is invalidated each time a project anywhere on the system is updated by a signal in cbh_core_model
        """
        if hasattr(request, "session"):
            serialized = request.session.get("projects_list_cache", None)
            if serialized and (request.GET.get("do_cache", False) or kwargs.get("do_cache", False)):
                desired_format = self.determine_format(request)
                print "from cache"
                return HttpResponse(content=serialized,
                                content_type=build_content_type(desired_format))
        return super(ChemregProjectResource, self).get_list(request, **kwargs)
        

    def post_list(self, request, **kwargs):
        """Ensure that when saving the project, any integrity errors from repeated names are passed up as a 409 conflict"""

        from django.db import IntegrityError
        try:
            return super(ChemregProjectResource, self).post_list(request, **kwargs)
        except IntegrityError, e:
            raise ImmediateHttpResponse(HttpConflict("Project with that name already exists"))


    def patch_detail(self, request, **kwargs):
        """Ensure that when saving the project, any integrity errors from repeated names are passed up as a 409 conflict"""

        from django.db import IntegrityError
        try:
            return super(ChemregProjectResource, self).patch_detail(request, **kwargs)
        except IntegrityError, e:
            raise ImmediateHttpResponse(HttpConflict("Project with that name already exists"))



    def patch_list(self, request, **kwargs):
        """Ensure that when saving the project, any integrity errors from repeated names are passed up as a 409 conflict"""

        from django.db import IntegrityError
        try:
            return super(ChemregProjectResource, self).patch_list(request, **kwargs)
        except IntegrityError, e:
            raise ImmediateHttpResponse(HttpConflict("Project with that name already exists"))




    def get_object_list(self, request):
        """Sort the projects by modified date descending by default"""
        return super(ChemregProjectResource,
                     self).get_object_list(request).prefetch_related(Prefetch('project_type'
                                                                              )).order_by('-modified')


    def alter_list_data_to_serialize(self, request, bundle):
        '''Here we append a list of tags to the data of the GET request if the
        search fields are required'''

        userres = UserResource()
        userbundle = userres.build_bundle(obj=request.user,
                                          request=request)
        userbundle = userres.full_dehydrate(userbundle)
        bundle['user'] = userbundle.data
        self._meta.authorization.alter_project_data_for_permissions(bundle, request)
        

        names_list = [] 
        #The tabular data schema is produced here as we have a list of the fields that the user has access to here
        schema = settings.TABULAR_DATA_SETTINGS["schema"].copy()

        #set up a restricted fields lookup for removing restricted fields from output
        bundle["user_restricted_fieldnames"] = {proj.data["id"] : proj.data["users_restricted_fields"] for proj in bundle[self._meta.collection_name]}

        for projbund in bundle["objects"]:
            for bund in projbund.data["custom_field_config"].data["project_data_fields"]:
                if bund.data["handsontable"]["data"] not in schema:
                    #The handsontable part is the same for each project
                    schema[bund.data["handsontable"]["data"]] = bund.data["handsontable"].copy()
                    names_list.append(bund.data["handsontable"]["data"])
                    schema[bund.data["handsontable"]["data"]]["project_specific_schema"] = {}
                    schema[bund.data["handsontable"]["data"]]["projects"] = []

                
                    #We have set up the app so there is one field per project so here we are just
                    # ensuring we know the field type, everything else can be looked up from 
                real_renderer_for_this_project = bund.obj.FIELD_TYPE_CHOICES[bund.obj.field_type]["data"]["renderer_named"]
                schema[bund.data["handsontable"]["data"]]["project_specific_schema"][projbund.data["id"]] = { "field_type" : bund.obj.field_type, "renderer_named" : real_renderer_for_this_project , "open_or_restricted" : bund.obj.open_or_restricted}

                del bund.data["handsontable"]
                #del bund.data["users_restricted_fields"]

        bundle["tabular_data_schema"] = {}
        bundle["tabular_data_schema"]["included_in_tables"] = {}
        for key, value in settings.TABULAR_DATA_SETTINGS.items():
            if key != "schema":
                bundle["tabular_data_schema"]["included_in_tables"][key] = {"default" : value["start"] + names_list + value["end"]}

        bundle["tabular_data_schema"]["schema"] = schema

        return bundle


    def alter_detail_data_to_serialize(self, request, bundle, tabular_data_schema=False):
        """
        Currently does not do all the things done by alter list data and perhaps better not to
        as this function is used to update a project 
        """
        userres = UserResource()
        userbundle = userres.build_bundle(obj=request.user,
                                          request=request)
        userbundle = userres.full_dehydrate(userbundle)
        bundle.data['user'] = userbundle.data

        if request.GET.get("tabular_data_schema", False):
            data = {"objects" : [bundle,]}
            self.alter_list_data_to_serialize( request, data)
            bundle.data = data["objects"][0].data
            bundle.data["tabular_data_schema"] = data["tabular_data_schema"]
        else:
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
        Ensures that the cache is invalidated and re saved as appropriate
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
        elif not hasattr(data, "data"):
            if data.get("objects", False) and (request.GET.get("do_cache", False) or response_kwargs.get("do_cache", False)) and response_class == HttpResponse:
                if hasattr(request, "session"):
                    request.session["projects_list_cache"] = serialized
        if response_class in [http.HttpCreated, http.HttpAccepted]:
            if hasattr(request, "session"):
                if "projects_list_cache" in request.session:
                    del request.session["projects_list_cache"]
            remove_session_cached_projectlists()

            create_or_update_project_index(data.obj.id)
        return rc
    
def create_or_update_project_index(project_id):
    project_json_with_tabular_schema = get_schemata([project_id,])
    for pd in project_json_with_tabular_schema:
        #There is only one pd!
        pass
    index_name = elasticsearch_client.get_project_index_name(project_id)
    elasticsearch_client.create_or_update_index(index_name, pd["tabular_data_schema"]["for_indexing"])


def get_schemata(project_ids, fieldlist_name="indexing", request=None):
    if project_ids is None:
        project_ids = get_model("cbh_core_model","Project").objects.filter().values_list("id", flat=True)
    crp = ChemregProjectResource()
    request_factory = RequestFactory()
    user = get_user_model().objects.filter(is_superuser=True)[0]
    
 
    for pid in project_ids:
        req = request_factory.get("/")
        if request:

            #Ensures that the restricted fields are obeyed if downloading data
            #Note that we cannot use the request raw because it will try to export XLSX because
            #of the way that tastypie works
            req.user = request.user
        else:
            #If this request is coming from the reindex_all_compounds request then we need to
            #Ensure that a superuser is aganst the request so we have access to all fields
            req.user = user
        req.GET = req.GET.copy()
        req.GET["tabular_data_schema"] = True

        data = crp.get_detail(req, pk=pid)

        projdata = json.loads(data.content)
        projdata["tabular_data_schema"]["for_%s" % fieldlist_name] = [ projdata["tabular_data_schema"]["schema"][field] 
                                  for field in projdata["tabular_data_schema"]["included_in_tables"][fieldlist_name]["default"]]
        yield projdata


def get_indexing_schemata(project_ids, fieldlist_name="indexing"):
    """Get a cached version of the project schemata in the right format 
    for elasticsearch indexing or for retrieving the Excel backup of the data"""

    proj_datasets = get_schemata(project_ids, fieldlist_name="indexing")
        #Do the dict lookups now to avoid multiple times later
    projdatadict = {}
    for projdata in proj_datasets:
        projdatadict[projdata["resource_uri"]] = projdata

    return projdatadict

class NoCustomFieldsChemregProjectResource(ChemregProjectResource):
    """A class for use in related resources in cases where 
    we do not want to save all or retrieve all of the custom fields such as when 
    indexing the CBHCompoundBatch objects"""
    custom_field_config = fields.ForeignKey(ChemRegCustomFieldConfigResource,
                                            'custom_field_config', blank=False, null=False)
    class Meta(ChemregProjectResource.Meta):
        resource_name = 'cbh_projects_nocfc'




class CBHNoSerializedDictField(fields.ApiField):
    """
    A dictionary field.
    Convert is no longer needed, we cast the data into the elasticsearch format elsewhere
    """
    dehydrated_type = 'dict'
    help_text = "A dictionary of data. Ex: {'price': 26.73, 'name': 'Daniel'}"




class CSRFExemptMixin(object):
    """CSRF exception to ensure that the user can press the back button after logging in
    and get the expected result"""
    @method_decorator(csrf_exempt)
    def dispatch(self, *args, **kwargs):
        return super(CSRFExemptMixin, self).dispatch(*args, **kwargs)

def get_field_name_from_key(key):
    """Possibly deprecated function to get a new field name when indexing data"""
    return key.replace(u"__space__", u" ")


def get_key_from_field_name(name):
    """Possibly deprecated function to get a new field name when indexing data"""
    return name.replace(u" ", u"__space__")




class SimpleResourceURIField(fields.ApiField):

    """
    Provide just the id field as a resource URI
    """
    dehydrated_type = 'string'
    is_related = False
    self_referential = False
    help_text = 'A related resource. Can be either a URI or set of nested resource data.'

    def __init__(self, to, attribute, full=False, related_name=None, default=fields.NOT_PROVIDED, null=False, blank=False, readonly=False,  unique=False, help_text=None, use_in='all'):
        """
        Init method very similar to a realtedresource in tastypie
        """
        super(SimpleResourceURIField, self).__init__(attribute=attribute, default=default, null=null, blank=blank,
                                                     readonly=readonly, unique=unique, help_text=help_text, use_in=use_in)
        self.related_name = related_name
        self.to = to
        self._to_class = None
        self._rel_resources = {}

        self.api_name = None
        self.resource_name = None

        if self.to == 'self':
            self.self_referential = True



    def contribute_to_class(self, cls, name):
        """
        method very similar to a realtedresource in tastypie
        """
        super(SimpleResourceURIField, self).contribute_to_class(cls, name)

        # Check if we're self-referential and hook it up.
        # We can't do this quite like Django because there's no ``AppCache``
        # here (which I think we should avoid as long as possible).
        if self.self_referential or self.to == 'self':
            self._to_class = cls

    def convert(self, value):
        """
        Handles conversion between the data found and the type of the field.
        Extending classes should override this method and provide correct
        data coercion.
        """
        if value is None:
            return None
        cls = self.to_class()
        resource_uri = cls.get_resource_uri()
        return "%s/%d" % (resource_uri, value)

    def hydrate(self, bundle):
        """
        Takes data stored in the bundle for the field and returns it. Used for
        taking simple data and building a instance object.
        """
        if self.readonly:
            return None
        if self.instance_name not in bundle.data:
            if self.is_related and not self.is_m2m:
                # We've got an FK (or alike field) & a possible parent object.
                # Check for it.
                if bundle.related_obj and bundle.related_name in (self.attribute, self.instance_name):
                    return bundle.related_obj
            if self.blank:
                return None
            if self.attribute:
                try:
                    val = getattr(bundle.obj, self.attribute, None)

                    if val is not None:
                        return val
                except ObjectDoesNotExist:
                    pass
            if self.instance_name:
                try:
                    if hasattr(bundle.obj, self.instance_name):
                        return getattr(bundle.obj, self.instance_name)
                except ObjectDoesNotExist:
                    pass
            if self.has_default():
                if callable(self._default):
                    return self._default()

                return self._default
            if self.null:
                return None

            raise ApiFieldError(
                "The '%s' field has no data and doesn't allow a default or null value." % self.instance_name)
        # New code to rerturn URI
        value = bundle.data[self.instance_name]
        if value is None:
            return value
        if str(value).endswith("/"):
            value = value[:-1]
        data = str(value).split("/")
        
        return int(data[len(data) - 1])
        
    @property
    def to_class(self):
        """
         method very similar to a realtedresource in tastypie
        """
        # We need to be lazy here, because when the metaclass constructs the
        # Resources, other classes may not exist yet.
        # That said, memoize this so we never have to relookup/reimport.
        if self._to_class:
            return self._to_class

        if not isinstance(self.to, six.string_types):
            self._to_class = self.to
            return self._to_class

        # It's a string. Let's figure it out.
        if '.' in self.to:
            # Try to import.
            module_bits = self.to.split('.')
            module_path, class_name = '.'.join(
                module_bits[:-1]), module_bits[-1]
            module = importlib.import_module(module_path)
        else:
            # We've got a bare class name here, which won't work (No AppCache
            # to rely on). Try to throw a useful error.
            raise ImportError(
                "Tastypie requires a Python-style path (<module.module.Class>) to lazy load related resources. Only given '%s'." % self.to)

        self._to_class = getattr(module, class_name, None)

        if self._to_class is None:
            raise ImportError("Module '%s' does not appear to have a class called '%s'." % (
                module_path, class_name))

        return self._to_class




class Index(TemplateView):
    """
    Serves the static template from the angularjs folder (the dist folder is the location that grunt build puts the production artifact
    We then copy this dist folder across with collectstatic by declaring it in the list of static folders
    """

    template_name = 'dist/index.html'  # or define get_template_names()

    def get(self, request, *args, **kwargs):
        """Ensure that the csrf token is always set in the cookie when rendering"""
        context = self.get_context_data(**kwargs)
        from django.middleware.csrf import get_token
        csrf_token = get_token(request)
        return self.render_to_response(context)



#-------------------------------------------------------------------------

class ProjectPermissionResource(ModelResource):
    """
    Allows updating of user permissions - project owners can change who views, edits and owns their project
    Data is retrieved using the codename of the 
    """
    users = fields.ToManyField("cbh_core_api.resources.UserResource",attribute="user_set")
    # groups = fields.ToManyField(GroupResource, attribute="group_set")
    codename = fields.CharField(readonly=True)

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<codename>[\w\d_.-]+)/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]

    class Meta:
    
        queryset = Permission.objects.all()
        resource_name = 'cbh_permissions'
        allowed_methods = ["get", "post", "patch", "put"]
        authentication = SessionAuthentication()
        authorization = ProjectPermissionAuthorization()
        detail_uri_name = 'codename'



    def save(self, bundle, skip_errors=False):
        """Ensure that the m2m bundle is hydrated before continuing as this is needed for the authorization"""
        if bundle.via_uri:
            return bundle

        self.is_valid(bundle)

        if bundle.errors and not skip_errors:
            raise ImmediateHttpResponse(response=self.error_response(bundle.request, bundle.errors))
        m2m_bundle = self.hydrate_m2m(bundle)

        # Check if they're authorized.
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), m2m_bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), m2m_bundle)

        # Save FKs just in case.
        self.save_related(bundle)

        # Save the main object.
        obj_id = self.create_identifier(bundle.obj)

        if obj_id not in bundle.objects_saved or bundle.obj._state.adding:
            bundle.obj.save()
            bundle.objects_saved.add(obj_id)

        # Now pick up the M2M bits.
        self.save_m2m(m2m_bundle)
        return bundle


#-------------------------------------------------------------------------



class Login( CSRFExemptMixin, FormView):
    """
    Login FormView using the standard locations for login templates in django
    """
    form_class = AuthenticationForm
    template_name = "cbh_chem_api/login.html"
    logout = None

    def get_context_data(self, **kwargs):
        context = super(Login, self).get_context_data(**kwargs)
        context["skinningconfig"] = SkinningConfig.objects.all()[0]
        return context
	
    def get(self, request, *args, **kwargs):
        """
        Ensure that the webauth is setup if required 
        Set the CSRF cookie
        """
        from django.middleware.csrf import get_token
        csrf_token = get_token(request)
        context = self.get_context_data(
            form=self.get_form(self.get_form_class()))
        redirect_to = settings.LOGIN_REDIRECT_URL

        if "django_webauth" in settings.INSTALLED_APPS:
            context["webauth_login"] = True
            username = request.META.get('REMOTE_USER', None)
            if not username:
                # Here we check if this was a redirect after logout in which
                # case we show the button to log out of webauth entirely
                username = request.META.get('HTTP_X_WEBAUTH_USER', None)
            if username:
                context["logout"] = True
        else:
            context["password_login"] = True

        if request.user.is_authenticated():
            return HttpResponseRedirect(redirect_to)
        return self.render_to_response(context)


    def post(self, request, *args, **kwargs):
        """
        Handles POST requests, instantiating a form instance with the passed
        POST variables and then checked for validity.
        """
        redirect_to = settings.LOGIN_REDIRECT_URL
        '''Borrowed from django base detail view'''
        from django.middleware.csrf import get_token
        csrf_token = get_token(request)
        if request.user.is_authenticated():
            #The user has pressed back in their browser and therefore should be redirected
            return HttpResponseRedirect(redirect_to)
        form = self.get_form(self.get_form_class())
        if form.is_valid():
            return self.form_valid(form)
        else:
            return self.form_invalid(form)



    def form_valid(self, form):
        """Remove test cookies as required"""
        redirect_to = settings.LOGIN_REDIRECT_URL

        auth_login(self.request, form.get_user())
        if self.request.session.test_cookie_worked():
            self.request.session.delete_test_cookie()
        # return self.render_to_response(self.get_context_data())
        return HttpResponseRedirect(redirect_to)

    def form_invalid(self, form):
        """Send the form back again if invalid with errors"""
        return self.render_to_response(self.get_context_data(form=form))



class Logout(View):
    """Standard logout view to just redirect to the login page"""
    def get(self, request, *args, **kwargs):
        auth_logout(request)
        return HttpResponseRedirect(settings.LOGOUT_REDIRECT_URL)


def build_content_type(format, encoding='utf-8'):
    """
    Appends character encoding to the provided format if not already present.
    """
    if 'charset' in format:
        return format

    return "%s; charset=%s" % (format, encoding)




class SkinningResource(ModelResource):

    '''URL resourcing for pulling out sitewide skinning config 
    Including the schemat for use on the search page
    '''
    tabular_data_schema = fields.DictField(default={}, help_text="empty dictionary, filled on the client side with the tabular data information in order to render handsontable")
    query_schemaform = fields.DictField(help_text="The angular schema form schema and form used to render the search forms in ChemBio Hub Platform")
    chem_query_schemaform = fields.DictField(help_text="The chemical angular schema form for chemical search")
    sort_schemaform = fields.DictField(help_text="The angular schema form json schema used for sorting fields")
    hide_schemaform = fields.DictField(help_text="The angular schema form json schema used for hiding fields")
    savedsearch_schemaform = fields.DictField(help_text="The angular schema form json schema used for the saved search feature")
    filters_applied = fields.ListField(default=[], help_text="empty list which is filled on the front end with the filters currently active on the dataset")
    sorts_applied = fields.ListField(default=[], help_text="empty list which is filled on the front end with the sorts currently active on the dataset")
    hides_applied = fields.ListField(default=[], help_text="empty list which is filled on the front end with the hidden fields on the dataset")

    filters_objects = fields.ListField(default=[], help_text="empty list which is filled on the front end with the columns which are currently filtered")
    sort_objects = fields.ListField(default=[], help_text="empty list which is filled on the front end with the columns which are currently sorted")
    hide_objects = fields.ListField(default=[],  help_text="empty list which is filled on the front end with the columns which are currently hidden")
    field_type_choices = fields.ListField(default=[],  help_text="List of the fields types which the user must pick from when creating a project")

    class Meta:
        always_return_data = True
        queryset = SkinningConfig.objects.all()
        resource_name = 'cbh_skinning'
        #authorization = Authorization()
        include_resource_uri = True
        allowed_methods = ['get', 'post', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()


    def dehydrate_savedsearch_schemaform(self, bundle):
        """Pull out the angular schema form json schema which is used when saving a search from the settings"""
        return settings.SAVED_SEARCH_SCHEMAFORM

    def dehydrate_query_schemaform(self, bundle):
        """Pull out the search form angular schema form json schema from the settings"""
        return settings.CBH_QUERY_SCHEMAFORM

    def dehydrate_hide_schemaform(self, bundle):
        """Pull out the angular schema form json schema from the settings which is used to hide a field"""
        return settings.CBH_HIDE_SCHEMAFORM

    def dehydrate_sort_schemaform(self, bundle):
        """Pull out the angular schema form json schema from the settings which is used to sort a field"""
        return settings.CBH_SORT_SCHEMAFORM

    def dehydrate_chem_query_schemaform(self, bundle):
        """Pull out the angular schema form json schema from the settings which is used to run chemical search"""
        return settings.CBH_CHEMICAL_QUERY_SCHEMAFORM

    def dehydrate_field_type_choices(self, bundle):
        """Pull out the field type choices from the models file which are used when adding a project
        todo consider moving project schemaform to the back end"""
        return [{"name": value["name"], "value": key} for key, value in PinnedCustomField.FIELD_TYPE_CHOICES.items()]



class TemplateProjectFieldResource(ModelResource):

    """Provides the schema information about a field that is required by front end apps"""
 
    class Meta:
        queryset = get_model("cbh_core_model","PinnedCustomField").objects.all()
        always_return_data = True
        resource_name = 'cbh_template_fields'
        include_resource_uri = False
        allowed_methods = ['get']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        authorization = Authorization()

    def dehydrate_id(self, bundle):
        return None



def get_field_list(project_type_bundle):
    """Given a project type object bundle return a project template 
    The project template is then used in the project adding form when the user clicks on that project type.
    The default project template is a single empty text field"""
    if project_type_bundle.obj.saved_search_project_type:
        return project_type_bundle.obj.SAVED_SEARCH_TEMPLATE
    elif project_type_bundle.obj.plate_map_project_type:
        return project_type_bundle.obj.PLATE_MAP_TEMPLATE
    elif len(project_type_bundle.data["custom_field_config_template"]) > 0:
        for field in  project_type_bundle.data["custom_field_config_template"]:
            field.data["id"] = None
        return [field.data for field in project_type_bundle.data["custom_field_config_template"]]
    else:
        return project_type_bundle.obj.DEFAULT_TEMPLATE

def get_fields(bundle):
    """May be a merge error need to lookj into this todo"""
    if bundle.obj.custom_field_config_template_id is None:
        return []
    return get_model("cbh_core_model","PinnedCustomField").objects.filter(custom_field_config_id=bundle.obj.custom_field_config_template_id)


class ProjectTypeResource(ModelResource):
    '''Resource for Project Type, specifies whether this is a chemical/inventory instance etc '''
    copy_action_name = fields.CharField(default="Clone", 
                                        help_text="The name for how the user clones the object changes dependent on the project type, this name comes from the back end and is used in front end buttons")
    custom_field_config_template = fields.ToManyField("cbh_core_api.resources.TemplateProjectFieldResource", 
        attribute=lambda bundle: get_fields(bundle) ,  
        full=True, 
        readonly=True, 
        null=True,
        help_text="Template set of fields for the projects that use this project type. Acts like a class inheritance")
    project_template = fields.DictField(default={}, help_text="Full project template")

    def alter_list_data_to_serialize(self, request, data):
        """Based on the custom field config template, generate a project template that can be used in the front end"""
        for bun in data["objects"]:
            bun.data["project_template"] =  {
                "project_type": bun.data["resource_uri"],
                    "custom_field_config": {
                        "project_data_fields": get_field_list(bun),
                        "name": ""
                    },
                    "name": ""
                }
        return data

    def dehydrate_copy_action_name(self, bundle):
        """Set the copy action name for the buttons based on the project type"""
        if bundle.obj.show_compounds:
            return "Clone /<br/> Add Structure"
        else:
            return "Clone Item"



    class Meta:
        """Setup for the project type resource. Be aware that anyone can currently update project types but we do not enable this feature in the front end"""
        always_return_data = True
        queryset = ProjectType.objects.all()
        resource_name = 'cbh_project_types'
        authorization = Authorization()
        include_resource_uri = True
        allowed_methods = ['get', 'post', 'patch', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        filtering = {
            "saved_search_project_type": ALL,
            "plate_map_project_type": ALL,
            "name": ALL
        }



class MyPasswordResetForm(PasswordResetForm):
    """Password reset form used both for password resets and for users who have just been invited and need to set a password"""
    def save(self, domain_override=None,
             subject_template_name='registration/password_reset_subject.txt',
             email_template_name='registration/password_reset_email.html',
             use_https=False, token_generator=default_token_generator,
             from_email=None, request=None, html_email_template_name=None, extra_email_context={}, user=None):
        """
        Generates a one-use only link for resetting password and sends to the
        user.
        """
        from django.core.mail import send_mail
        email = self.cleaned_data["email"]
 
        if not domain_override:
            #The sites framework does not necessarily give the correct URL of the deployed system so we get it from the request instead
            current_site = get_current_site(request)
            site_name = current_site.name
            domain = current_site.domain
        else:
            site_name = domain = domain_override
        c = {
            'email': user.email,
            'domain': domain,
            'site_name': site_name,
            'uid': urlsafe_base64_encode(force_bytes(user.pk)),
            'user': user,
            'token': token_generator.make_token(user),
            'protocol': 'https' if use_https else 'http',
            'extra' : extra_email_context,

        }
        subject = loader.render_to_string(subject_template_name, c)
        # Email subject *must not* contain newlines
        subject = ''.join(subject.splitlines())
        html_email = loader.render_to_string(email_template_name, c)
        soup = BeautifulSoup(html_email)
        email = soup.getText()

        
        send_mail(subject, email, from_email, [user.email], html_message=html_email, fail_silently=False)

 

class InvitationResource(UserHydrate, ModelResource):
    '''Resource for Invitation model. This will setup creation of the invite email and new user '''

    created_by = fields.ForeignKey(
        "cbh_core_api.resources.UserResource", 'created_by', full=True, help_text="creator of the invitation")
    class Meta:
        queryset = Invitation.objects.all()
        resource_name = 'invitations'
        authorization = InviteAuthorization()
        include_resource_uri = True
        allowed_methods = ['get', 'post', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        always_return_data = True
        filtering = {
            "email": ALL_WITH_RELATIONS
        }





    def get_form(self, email, new_user, data, created, request, email_template_name, subject_template_name):
        """Save the data from the invitation form (results in an email being sent)"""
        server = settings.SERVER_EMAIL
        form = MyPasswordResetForm(QueryDict(urlencode({"email": email})))
        hostname = request.META["HTTP_ORIGIN"]
        if form.is_valid():
            form.users_cache = [new_user,]
            opts = {
                'use_https': request.is_secure(),
                'token_generator': default_token_generator,
                'from_email': server,
                'user' : new_user,
                'email_template_name': email_template_name,
                'subject_template_name': subject_template_name,
                'request': request,
                'extra_email_context': {'hostname':hostname, 'invite': data.data, 'login_url' : settings.LOGIN_URL, },
            }
            form.save(**opts)

        else:
            raise BadRequest("Email not valid")

    def create_response(self, request, data, response_class=HttpResponse, **response_kwargs):
        """
        Extracts the common "which-format/serialize/return-response" cycle.
        Mostly a useful shortcut/hook.
        """
        desired_format = self.determine_format(request)
       

        if response_class == http.HttpCreated:     
            email = data.data["email"]
            if email.endswith("ox.ac.uk"):
                #send via webauth
                raise BadRequest("We do not yet support inviting users at Oxford to projects. This feature will come soon.")
            else:

                UserObj = get_user_model()
                new_user, created = UserObj.objects.get_or_create(email=email, username=email)
                logger.info(data.data)
                for perm in data.data["projects_selected"]:
                    p = Project.objects.get(id=perm["id"])
                    p.make_viewer(new_user)
                    p.save()
                data.data["message"] = "Invite sent successfully to %s, would you like to invite anyone else?" % email
                email_template_name = 'cbh_core_api/email_new_user.html'
                subject_template_name = 'cbh_core_api/subject_new_user.html'
                if not created:
                    projects_with_reader_access = viewer_projects(new_user)
                    all_projects_equal = True

                    all_selected_ids = set([new_proj["id"] for new_proj in data.data["projects_selected"]])
                    new_ids = all_selected_ids - set(projects_with_reader_access)
                    
                    if(len(new_ids) > 0):
                        email_template_name = 'cbh_core_api/email_project_access_changed.html'
                        subject_template_name = 'cbh_core_api/subject_project_access_changed.html'
                        all_projects_equal = False
                        data.data["message"] = "Existing user %s invited to new projects, would you like to invite anyone else?" % email
                    else:
                        if not data.data.get("remind", False):
                            raise ImmediateHttpResponse(http.HttpConflict('{"error": "User already exists, do you wish to invite again?"}'))
                        if new_user.has_usable_password():
                            email_template_name = 'cbh_core_api/email_reminder.html'
                            subject_template_name = 'cbh_core_api/subject_reminder.html'
                            data.data["message"] = "Sign-up reminder sent to %s, would you like to invite anyone else?" % email
                        else:
                            email_template_name = 'cbh_core_api/email_reminder_already_logged_on.html'
                            subject_template_name = 'cbh_core_api/subject_reminder.html'
                            data.data["message"] = "User %s reminded to look at these projects, would you like to invite anyone else?" % email
                form = self.get_form( email, new_user, data, created, request, email_template_name, subject_template_name)         

        serialized = self.serialize(request, data, desired_format)
        rc = response_class(content=serialized, content_type=build_content_type(
            desired_format), **response_kwargs)       
        return rc



class UserResource(ModelResource):
    '''Displays information about the User's privileges and personal data'''
    can_view_chemreg = fields.BooleanField(default=True, help_text="Whether the user has chemireg enabled deprecated")
    can_view_assayreg = fields.BooleanField(default=True, help_text="Whether the user has assayreg enabled deprecated")
    is_logged_in = fields.BooleanField(default=False,  help_text="Whether this user in the list is the logged in user")
    can_create_and_own_projects = fields.BooleanField(default=False,  help_text="Whether this user is allowed to create or own projects")
    display_name = fields.CharField(default="", help_text="How we want the name of this user to be displayed on the front end")

    class Meta:
        filtering = {
            "username": ALL_WITH_RELATIONS,
            "id": ALL_WITH_RELATIONS,
        }
        queryset = User.objects.all()
        resource_name = 'users'
        allowed_methods = ["get", "post"]
        excludes = ['email', 'password', 'is_active']
        authentication = SessionAuthentication()
        authorization = Authorization()


    def apply_authorization_limits(self, request, object_list):
        """probably deprecated"""
        return object_list.get(pk=request.user.id)

    def get_object_list(self, request):
        """
        deprecated override
        """
        return super(UserResource, self).get_object_list(request)


    def dehydrate_display_name(self, bundle):
        """How we want the username or first and last name to be displayed on the system"""
        if bundle.obj.first_name:
            return "%s %s" % (bundle.obj.first_name, bundle.obj.last_name)
        else:
            return bundle.obj.username

    def dehydrate_can_create_and_own_projects(self, bundle):
        """Internal users (denoted by their email pattern match) are allowed to add and own projects"""
        if bundle.obj.is_superuser:
            return True
        perms = bundle.obj.get_all_permissions()
        if "cbh_core_model.add_project" in perms:
            return True
        return False


    def dehydrate_is_logged_in(self, bundle):
        """
        Whether the user is logged in
        """
        if hasattr(bundle.request, "user"):
            if bundle.obj.id == bundle.request.user.id:
                return True
        return False

    def dehydrate_can_view_chemreg(self, bundle):
        '''The cbh_core_model.no_chemreg role in the Django admin is used to
        deny access to chemreg. As superusers have all permissions  by 
        default they would be denied access therefore we check for superuser status and allow access'''
        if bundle.obj.is_superuser:
            return True
        perms = bundle.obj.get_all_permissions()
        if "cbh_core_model.no_chemreg" in perms:
            return False
        return True

    def dehydrate_can_view_assayreg(self, bundle):
        '''The cbh_core_model.no_assayreg role in the Django admin is used to
        deny access to assayreg. As superusers have all permissions  by 
        default they would be denied access therefore we check for superuser status and allow access'''

        if bundle.obj.is_superuser:
            return True
        perms = bundle.obj.get_all_permissions()
        if "cbh_core_model.no_assayreg" in perms:
            return False
        return True



def chunk_upload_to(instance, filename):
    """
    Save chunk to the right path and filename based in is number
    """
    return os.path.join(FLOWJS_PATH, instance.filename)


def remove_expired_files():
    """
    Remove non completed uploads
    """
    from cbh_core_model.models import FlowFile
    FlowFile.objects.filter(
        state__in=[FlowFile.STATE_UPLOADING, FlowFile.STATE_UPLOAD_ERROR],
        updated__lte=datetime.datetime.date() - datetime.timedelta(days=FLOWJS_EXPIRATION_DAYS)
    ).delete()
