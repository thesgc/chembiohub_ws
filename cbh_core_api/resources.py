import logging

# Get an instance of a logger
logger = logging.getLogger(__name__)


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
from cbh_core_model.models import DataType
from cbh_core_model.models import Project
from cbh_core_model.models import ProjectType
from cbh_core_model.models import SkinningConfig
from cbh_core_model.models import Invitation


from cbh_core_api.authorization import ProjectListAuthorization, InviteAuthorization, viewer_projects, ProjectPermissionAuthorization
from tastypie.authentication import SessionAuthentication
from tastypie.paginator import Paginator
from cbh_core_api.serializers import CustomFieldsSerializer

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

class CBHNoSerializedDictField(fields.ApiField):
    """
    A dictionary field.
    """
    dehydrated_type = 'dict'
    help_text = "A dictionary of data. Ex: {'price': 26.73, 'name': 'Daniel'}"

    def convert(self, value):
        data = dict(value)
        for key, v in data.iteritems():
            if isinstance(v, basestring):
                if v.startswith("[") and v.endswith("]"):
                    try:
                        data[k] = json.loads(v)
                        continue
                    except:
                        pass
                if v.startswith("{") and v.endswith("}"):
                    try:
                        data[k] = json.loads(v)
                        continue
                    except:
                        pass
        return data




class UserHydrate(object):
    def hydrate_created_by(self, bundle):
        User = get_user_model()

        if bundle.obj.id:
            pass
        else:
            user = User.objects.get(pk=bundle.request.user.pk)
            bundle.obj.created_by = user
        return bundle


class CSRFExemptMixin(object):
    @method_decorator(csrf_exempt)
    def dispatch(self, *args, **kwargs):
        return super(CSRFExemptMixin, self).dispatch(*args, **kwargs)

def get_field_name_from_key(key):
    return key.replace(u"__space__", u" ")


def get_key_from_field_name(name):
    return name.replace(u" ", u"__space__")


from django.middleware.csrf import get_token


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

    template_name = 'dist/index.html'  # or define get_template_names()

    def get(self, request, *args, **kwargs):
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
    form_class = AuthenticationForm
    template_name = "cbh_chem_api/login.html"
    logout = None

    def get(self, request, *args, **kwargs):

        from django.middleware.csrf import get_token
        csrf_token = get_token(request)
        context = self.get_context_data(
            form=self.get_form(self.get_form_class()))
        redirect_to = settings.LOGIN_REDIRECT_URL
        '''Borrowed from django base detail view'''
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
        redirect_to = settings.LOGIN_REDIRECT_URL

        auth_login(self.request, form.get_user())
        if self.request.session.test_cookie_worked():
            self.request.session.delete_test_cookie()
        # return self.render_to_response(self.get_context_data())
        return HttpResponseRedirect(redirect_to)

    def form_invalid(self, form):
        return self.render_to_response(self.get_context_data(form=form))

    # def dispatch(self, request, *args, **kwargs):
    #     request.session.set_test_cookie()
    #     return super(Login, self).dispatch(request, *args, **kwargs)


class Logout(View):

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

    '''URL resourcing for pulling out sitewide skinning config '''
    tabular_data_settings = fields.DictField()
    query_schema = fields.DictField()
    query_form = fields.DictField()

    class Meta:
        always_return_data = True
        queryset = SkinningConfig.objects.all()
        resource_name = 'cbh_skinning'
        #authorization = Authorization()
        include_resource_uri = True
        allowed_methods = ['get', 'post', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()

    def dehydrate_tabular_data_settings(self, bundle):
        return settings.TABULAR_DATA_SETTINGS

    def dehydrate_query_schema(self, bundle):
        return settings.CBH_QUERY_SCHEMA

    def dehydrate_query_form(self, bundle):
        return settings.CBH_QUERY_FORM



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
    if project_type_bundle.obj.saved_search_project_type:
        return project_type_bundle.obj.SAVED_SEARCH_TEMPLATE
    else:
        for field in  project_type_bundle.data["custom_field_config_template"]:
            field.data["id"] = None
        return [field.data for field in project_type_bundle.data["custom_field_config_template"]]



class ProjectTypeResource(ModelResource):
    '''Resource for Project Type, specifies whether this is a chemical/inventory instance etc '''
    copy_action_name = fields.CharField(default="Clone")
    custom_field_config_template = fields.ToManyField("cbh_core_api.resources.TemplateProjectFieldResource", attribute=lambda bundle: get_model("cbh_core_model","PinnedCustomField").objects.filter(custom_field_config_id=bundle.obj.custom_field_config_template_id) ,  full=True, readonly=True, null=True)
    project_template = fields.DictField(default={})

    def alter_list_data_to_serialize(self, request, data):
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
        if bundle.obj.show_compounds:
            return "Clone / Add Structure"
        else:
            return "Clone Item"



    class Meta:
        always_return_data = True
        queryset = ProjectType.objects.all()
        resource_name = 'cbh_project_types'
        authorization = Authorization()
        include_resource_uri = True
        allowed_methods = ['get', 'post', 'patch', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        filtering = {
            "saved_search_project_type": ALL
        }

class CustomFieldConfigResource(ModelResource):

    '''Resource for Custom Field Config '''
    class Meta:
        always_return_data = True
        queryset = CustomFieldConfig.objects.all()
        resource_name = 'cbh_custom_field_configs'
        #authorization = ProjectListAuthorization()
        include_resource_uri = True
        allowed_methods = ['get', 'post', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        filtering = {
            "name": ALL_WITH_RELATIONS
        }


class DataTypeResource(ModelResource):

    '''Resource for data types'''
    plural = fields.CharField(null=True)

    class Meta:
        always_return_data = True
        queryset = DataType.objects.all()
        resource_name = 'cbh_data_types'
        authorization = Authorization()
        include_resource_uri = True
        allowed_methods = ['get', 'post', 'patch', 'put']
        default_format = 'application/json'
        authentication = SessionAuthentication()
        filtering = {
            "name": ALL_WITH_RELATIONS
        }
        authorization = Authorization()

    def dehydrate_plural(self, bundle):
        return inflection.pluralize(bundle.obj.name)



class MyPasswordResetForm(PasswordResetForm):
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
        "cbh_core_api.resources.UserResource", 'created_by', full=True)
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


from django.contrib.auth.models import User

class UserResource(ModelResource):
    '''Displays information about the User's privileges and personal data'''
    can_view_chemreg = fields.BooleanField(default=True)
    can_view_assayreg = fields.BooleanField(default=True)
    is_logged_in = fields.BooleanField(default=False)
    can_create_and_own_projects = fields.BooleanField(default=False)
    display_name = fields.CharField(default="")

    class Meta:
        filtering = {
            "username": ALL_WITH_RELATIONS
        }
        queryset = User.objects.all()
        resource_name = 'users'
        allowed_methods = ["get", "post"]
        excludes = ['email', 'password', 'is_active']
        authentication = SessionAuthentication()
        authorization = Authorization()


    def apply_authorization_limits(self, request, object_list):
        return object_list.get(pk=request.user.id)

    def get_object_list(self, request):
        # return super(UserResource,
        # self).get_object_list(request).filter(pk=request.user.id)
        return super(UserResource, self).get_object_list(request)


    def dehydrate_display_name(self, bundle):
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

class CoreProjectResource(ModelResource):
    project_type = fields.ForeignKey(
        ProjectTypeResource, 'project_type', blank=False, null=False, full=True)
    custom_field_config = fields.ForeignKey(
        CustomFieldConfigResource, 'custom_field_config', blank=False, null=True, full=True)

    class Meta:
        queryset = Project.objects.all()
        authentication = SessionAuthentication()
        paginator_class = Paginator
        allowed_methods = ['get']
        resource_name = 'cbh_projects'
        authorization = ProjectListAuthorization()
        include_resource_uri = True
        default_format = 'application/json'
        #serializer = Serializer()
        serializer = CustomFieldsSerializer()
        filtering = {

            "project_key": ALL_WITH_RELATIONS,
        }

    def get_object_list(self, request):
        return super(CoreProjectResource, self).get_object_list(request).prefetch_related(Prefetch("project_type")).order_by('-modified')

    def alter_list_data_to_serialize(self, request, bundle):
        '''Here we append a list of tags to the data of the GET request if the
        search fields are required'''
        userres = UserResource()
        userbundle = userres.build_bundle(obj=request.user, request=request)
        userbundle = userres.full_dehydrate(userbundle)
        bundle['user'] = userbundle.data

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

#FlowFile relocation
import os
import datetime
from cbh_core_api.flowjs_settings import FLOWJS_PATH, FLOWJS_EXPIRATION_DAYS


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
