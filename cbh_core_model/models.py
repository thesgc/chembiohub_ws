# -*- coding: utf-8 -*-
"""Core models for ChemBio Hub platform, covering objects required for configuration of the tool such as projects and skinning"""
from django.db import models, connection

from solo.models import SingletonModel
from django_extensions.db.models import TimeStampedModel
from django.db.models.signals import post_save, pre_save
from django.contrib.contenttypes.models import ContentType
from django.contrib.auth.models import Permission, User, Group
from collections import OrderedDict
from django.utils.functional import cached_property
from copy import copy, deepcopy
import json
import dateutil
import time
import django
#from cbh_core_model.models import FlowFile
from django.template.defaultfilters import slugify 
from django_hstore import hstore
#FlowFile relocation stuff

import os
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from django.core.files.storage import default_storage
from django.conf import settings
from cbh_core_api.flowjs_settings import FLOWJS_PATH, FLOWJS_REMOVE_FILES_ON_DELETE, FLOWJS_AUTO_DELETE_CHUNKS
from cbh_core_api.utils import chunk_upload_to
from cbh_utils.idgenerator import IncrementalIdGenerator


PERMISSION_CODENAME_SEPARATOR = "__"
OPEN = "open"
RESTRICTED = "restricted"

RESTRICTION_CHOICES = ((OPEN, "Open to all viewers"),
                        (RESTRICTED, "Restricted to editors"))
#Assume that these are declared in order of increasing permission
PROJECT_PERMISSIONS = (("viewer", "Can View", {"linked_field_permission": RESTRICTED}),
                       ("editor", "Can edit or add batches",  {"linked_field_permission": OPEN}),
                       ("owner", "Can assign permissions",  {"linked_field_permission": OPEN}))



def get_all_project_ids_for_user_perms(perms, possible_perm_levels):
    """New implementation of permissions by project - the contenttype is left now as project and 
    the permission codename contains all information about the permission"""
    pids = []
    for perm in perms:
        proj_role = str(perm).split(".")[1]
        if PERMISSION_CODENAME_SEPARATOR in proj_role:
            prms = proj_role.split(PERMISSION_CODENAME_SEPARATOR)
            pid = prms[0]
            if pid[0].isdigit() and prms[1] in possible_perm_levels:
                pids.append(int(pid))
    return list(set(pids))


def get_old_project_ids_for_user_perms(perms, possible_perm_levels):
    """possibly deprecated"""
    pids = []
    for perm in perms:
        prms = str(perm).split(".")
        
        pid = prms[0]
        if pid[0].isdigit() and prms[1] in possible_perm_levels:
            pids.append(int(pid))
    return pids


def get_old_all_project_ids_for_user(user, possible_perm_levels):
    """deprecated"""
    return get_old_project_ids_for_user_perms(user.get_all_permissions(), possible_perm_levels)

def get_old_all_project_ids_for_user(user, possible_perm_levels):
    """deprecated"""
    return get_old_project_ids_for_user_perms(user.get_all_permissions(), possible_perm_levels)



def get_all_project_ids_for_user(user, possible_perm_levels):
    """Given a list of permission levels, extract all of the project ids for which the user has one of those permission levels"""
    return get_all_project_ids_for_user_perms(user.get_all_permissions(), possible_perm_levels)


def get_projects_where_fields_restricted(user):
    """Iterate through the permission choices available and assign a dictionary for this user of the index that should
    be used for each permission. 
    This is implemented in this way to allow the roles to be decopupled from the index that the user views"""
    indexes_dict = { OPEN : set(), RESTRICTED :set()}
    for perm in PROJECT_PERMISSIONS:
        for pid in get_all_project_ids_for_user(user, [perm[0]]):
            indexes_dict[perm[2]["linked_field_permission"]].add(pid)
    #Poen permission trump restricted ones
    indexes_dict[RESTRICTED] = indexes_dict[RESTRICTED] - indexes_dict[OPEN]
    return indexes_dict


class ProjectPermissionManager(models.Manager):
    """Manager methods to be inherited by the project object to manage permissions on projects"""
    def sync_all_permissions(self):
        """Run trough all of the projects on the system and sync up the permission objects that they require"""
        for proj in self.all():
            proj.sync_permissions()



    def get_next_incremental_id_for_compound(self, project_id):
        """Using the Redis ID generator retrieve the next ID for this project"""
        gen = IncrementalIdGenerator("project_%d" % project_id, maxReserveBuffer=1)
        new_id = gen.getId()
        if new_id == 1:
            proj = self.get(pk=project_id)
            while new_id != proj.project_counter_start:
                new_id = gen.getId()
        return new_id


def get_permission_name(name, permission):
    """Generate a permission name so that we have something to display in the django admin UI in the dropdown selector"""
    return "%s%s%s" % (name, PERMISSION_CODENAME_SEPARATOR, permission)

def get_permission_codename(id, permission):
    """For a given project a set of permissions are generated in django like
    2__owner
    2__editor
    2__viewer
    This function generates these codenames"""
    return "%d%s%s" % (id, PERMISSION_CODENAME_SEPARATOR, permission)


class ProjectPermissionMixin(models.Model):

    '''The aim of this mixin is to create a permission content type and a permission model for a given project
    It allows for pruning the contnet types once the model is changed
    '''

    objects = ProjectPermissionManager()


    def sync_permissions(self):
        '''first we delete the existing permissions that are not labelled in the model'''
        for perm in PROJECT_PERMISSIONS:
            self.get_instance_permission_by_codename(perm[0])




    def get_contenttype_for_instance(self):
        """All project permissions now simply sit under the project content type"""
        ct = ContentType.objects.get_for_model(self)
        return ct



    def get_instance_permission_by_codename(self, codename):
        """Get or create a permission object for a given codename on this project"""
        pm, created = Permission.objects.get_or_create(
            codename=get_permission_codename(self.id, codename), 
            content_type_id=self.get_contenttype_for_instance().id, 
            name=get_permission_name(self.name, codename))
        return pm

    def _add_instance_permissions_to_user_or_group(self, group_or_user, codename):
        """Give a user a given permission on a project"""
        if type(group_or_user) == Group:
            group_or_user.permissions.add(
                self.get_instance_permission_by_codename(codename))
        if type(group_or_user) == User:
            group_or_user.user_permissions.add(
                self.get_instance_permission_by_codename(codename))

    

    def make_editor(self, group_or_user):
        """Make a user object an editor of a project"""
        self._add_instance_permissions_to_user_or_group(
            group_or_user,  "editor")

    def make_viewer(self, group_or_user):
        """Make a user object an viewer of a project"""
        self._add_instance_permissions_to_user_or_group(
            group_or_user,   "viewer")

    def make_owner(self, group_or_user):
        """Make a user object an owner of a project"""
        self._add_instance_permissions_to_user_or_group(
            group_or_user, "owner")

    class Meta:
        """Abstract model as only used to inherit into project"""
        abstract = True


class ProjectType(TimeStampedModel):
    """The project type object"""
    plateSizes = "96,48"
    plateTypes = "working,backup"
    SAVED_SEARCH_TEMPLATE = [{
                            "required": False,
                            "field_type": "char",
                            "open_or_restricted": "open",
                            "name": "Alias"
                        },
                        {
                            "required": False,
                            "field_type": "char",
                            "open_or_restricted": "open",
                            "name": "URL"
                        }]
    PLATE_MAP_TEMPLATE = [
                        {"required":True, "field_type": "char", "open_or_restricted": "open", "name": "Name", 'readonly_on_create':True},
                        {"required":True, "field_type": "radios", "allowed_values": plateSizes, "open_or_restricted": "open", "name": "Plate Size",'readonly_on_create':True},
                        {"required":False, "field_type": "char", "open_or_restricted": "open", "name": "Description",'readonly_on_create':True},
                        {"required":True, "field_type": "radios", "allowed_values": plateTypes, "open_or_restricted": "open", "name": "Plate Type",'readonly_on_create':True},
                      ]
    DEFAULT_TEMPLATE = [{"required":False, "field_type": "char", "open_or_restricted": "open"},]
    ''' Allows configuration of parts of the app on a per project basis - initially will be used to separate out compound and inventory projects '''
    name = models.CharField(
        max_length=100, db_index=True, null=True, blank=True, default=None, help_text="Name of the project type")
    show_compounds = models.BooleanField(default=True, help_text="Whether this project type should show compounds")
    saved_search_project_type = models.BooleanField(default=False, help_text="Whether this project type is to be used in saved searches")
    plate_map_project_type = models.BooleanField(default=False, help_text="Whether this project type is for plate maps")
    #Set the default to -1 to match the empty default custom field config
    custom_field_config_template_id = models.IntegerField(default=None, null=True, blank=True, help_text="If you wish for this project type to include a set of default fields, find an appropriate custom field config ID to use as a template and add it here")
    set_as_default = models.BooleanField(default=False, help_text="Whether this project type should be used as the default project type on the system, set back to false if another project type is set as default")
    def __unicode__(self):
        return self.name


def make_default(sender, instance, **kwargs):
    """Ensure that there is only one default project type set at a time"""
    if instance.set_as_default:
        ProjectType.objects.exclude(pk=instance.id).update(set_as_default=False)


post_save.connect(make_default, sender=ProjectType, dispatch_uid="ptype")





class CustomFieldConfig(TimeStampedModel):
    '''
    Stores the schema of fields used in ChemiReg and AssayReg
    '''
    name = models.CharField(unique=True, max_length=500, null=False, blank=False)
    created_by = models.ForeignKey("auth.User")
    schemaform = models.TextField(default="", null=True, blank=True, )

    def __unicode__(self):
        return self.name

    def get_space_replaced_name(self):
        return self.name.replace(u" ", u"__space__")


class DataFormConfig(TimeStampedModel):

    '''deprecated Shared configuration object - all projects can see this and potentially use it
    Object name comes from a concatentaion of all of the levels of custom field config
    '''
    created_by = models.ForeignKey("auth.User")
    human_added = models.NullBooleanField(default=True)
    parent = models.ForeignKey(
        'self', related_name='children', default=None, null=True, blank=True)

    l0 = models.ForeignKey("cbh_core_model.CustomFieldConfig",
                           related_name="l0",
                           help_text="The first level in the hierarchy of the form you are trying to create. For example, if curating industries, companies,  employees , teams and departments, l0 would be industries.")
    l1 = models.ForeignKey("cbh_core_model.CustomFieldConfig",
                           related_name="l1",
                           null=True,
                           blank=True,
                           default=None,
                           help_text="The second level in the hierarchy of the form you are trying to create. For example, if curating industries, companies,  employees , teams and departments, l1 would be companies.")
    l2 = models.ForeignKey("cbh_core_model.CustomFieldConfig",
                           related_name="l2", null=True, blank=True, default=None,
                           help_text="The third level in the hierarchy of the form you are trying to create. For example, if curating industries, companies,  employees , teams and departments, l2 would be departments.")
    l3 = models.ForeignKey("cbh_core_model.CustomFieldConfig",
                           related_name="l3",
                           null=True,
                           blank=True,
                           default=None,
                           help_text="The forth level in the hierarchy of the form you are trying to create. For example, if curating industries, companies,  employees , teams and departments, l3 would be teams.")
    l4 = models.ForeignKey("cbh_core_model.CustomFieldConfig",
                           related_name="l4",
                           null=True,
                           blank=True,
                           default=None,
                           help_text="The fifth level in the hierarchy of the form you are trying to create. For example, if curating industries, companies,  employees , teams and departments, l4 would be employees.")

    def __unicode__(self):
        string = ""
        if self.l0:
            string += self.l0.__unicode__()
        if self.l1:
            string += " >> " + self.l1.__unicode__()
        if self.l2:
            string += " >> " + self.l2.__unicode__()
        if self.l3:
            string += " >> " + self.l3.__unicode__()
        if self.l4:
            string += " >> " + self.l4.__unicode__()
        return string

    class Meta:
        unique_together = (('l0', 'l1', 'l2', 'l3', 'l4'),)
        ordering = ('l0', 'l1', 'l2', 'l3', 'l4')

    def last_level(self):
        last_level = ""
        if self.l4_id is not None:
            return "l4"
        if self.l3_id is not None:
            return "l3"
        if self.l2_id is not None:
            return "l2"
        if self.l1_id is not None:
            return "l1"
        if self.l0_id is not None:
            return "l0"
        return last_level

    def get_all_ancestor_objects(obj, request, tree_builder={}, uri_stub=""):
        levels = ["l0", "l1", "l2", "l3", "l4"]
        levels = ["%s_id" % l for l in levels]
        used_levels = []
        for lev in levels:
            if getattr(obj, lev) is not None and lev != "%s_id" % obj.last_level():
                used_levels.append(lev)

        filters = []
        for i in range(1, len(used_levels)+1):
            level_name = used_levels[i-1]
            level_filters = {lev: getattr(obj, lev) for lev in used_levels[:i]}
            filters.insert(0, level_filters)

        for index, filter_set in enumerate(filters):
            defaults = {lev: None for lev in levels}
            new_filters = defaults.update(filter_set)

            defaults["defaults"] = {"created_by_id": request.user.id,
                                    "human_added": False}
            new_object, created = DataFormConfig.objects.get_or_create(
                **defaults)
            if not obj.parent_id:
                obj.parent_id = new_object.id
                obj.save()
#

            permitted_child_array = tree_builder.get(
                "%s/%d" % (uri_stub, new_object.id), [])
            permitted_child_array.append(obj)
            tree_builder[
                "%s/%d" % (uri_stub, new_object.id)] = list(set(permitted_child_array))

            obj = new_object
            if index == len(filters) - 1:
                tree_builder["root"] = [obj]


class Project(TimeStampedModel, ProjectPermissionMixin):
    ''' Project is a holder for moleculedictionary objects and for batches and permissions on the system are inherited from the project'''
    name = models.CharField(
        max_length=100, db_index=True, null=True, blank=True, default=None, help_text="What would you like to call your project")
    project_key = models.SlugField(
        max_length=50, db_index=True, null=True, blank=True, default=None, unique=True, help_text="key for the project, mostly deprecated in favour of id now but still used in some places")
    created_by = models.ForeignKey("auth.User", help_text="User who created the project")
    custom_field_config = models.ForeignKey(
        "cbh_core_model.CustomFieldConfig", related_name="project", null=True, blank=True, default=None, help_text="Custom field config object which is linked to this project")
    project_type = models.ForeignKey(
        ProjectType, null=True, blank=True, default=None, help_text="Project type that this project is linked to")
    project_counter_start = models.IntegerField(default=1, help_text="start of the incremental ID field for this project")

    class Meta:
        get_latest_by = 'created'

    def __unicode__(self):
        return self.name

    @models.permalink
    def get_absolute_url(self):
        return {'post_slug': self.project_key}


def get_name_for_custom_field_config_from_project(p):
    return "%d__%s__project__config" % (p.id, p.name)


def sync_permissions(sender, instance, created, **kwargs):
    '''After saving the project make sure it has entries in the permissions table
    After possible renaming of the project ensure that all associated Custom Field Configs etc. are renamed'''
    if created is True:
        instance.sync_permissions()

        instance.make_owner(instance.created_by)
    else:
        #Iterate all of the project permissions associated with this instance and update the user friendly name for those permissions
        proj_ct = ContentType.objects.get_for_model(instance)
        for perm_name in PROJECT_PERMISSIONS:
            perms = Permission.objects.filter(codename=get_permission_codename(instance.id, perm_name[0]), 
                content_type=proj_ct)
            for perm in perms:
                perm.name = get_permission_name(instance.name, perm_name[0])
                perm.save()

        instance.custom_field_config.name = get_name_for_custom_field_config_from_project(instance)
        instance.custom_field_config.save()

post_save.connect(sync_permissions, sender=Project, dispatch_uid="proj_perms")

def update_project_key(sender, instance, **kwargs):
    """Set the project key for a given project"""
    instance.project_key = slugify(instance.name)


pre_save.connect(update_project_key, sender=Project, dispatch_uid="proj_key")




class SkinningConfig(SingletonModel):

    '''Holds information about custom system messages and other customisable elements'''
    #created_by = models.ForeignKey("auth.User")
    instance_alias = models.CharField(
        max_length=50, null=True, blank=False, default='ChemiReg', help_text="What the instance of ChemBio Hub should be labelled as (not really used, deprecated")
    project_alias = models.CharField(
        max_length=50, null=True, blank=False, default='project', help_text="What projects on the system should be called (not really used, could be deprecated)")
    result_alias = models.CharField(
        max_length=50, null=True, blank=False, default='result', help_text="What search results on the system should be called (not really used, could be deprecated)")
    file_errors_from_backend = models.NullBooleanField(default=False, help_text="Whether the error message on file upload should be the exception message from the back end")
    enable_smiles_input = models.NullBooleanField(default=True, help_text="Whether the SMILES input field should be enabled when adding multiple compounds")
    data_manager_email = models.CharField(max_length=100, default="", help_text="Email address for support requests")
    data_manager_name = models.CharField(max_length=100, default="", help_text="Name of support requests manager")
    def __unicode__(self):
        return u"Skinning Configuration"

    class Meta:
        verbose_name = "Skinning Configuration"


def test_string(value):
    """Possibly deprecated"""
    return True

def test_bool(value):
    """Possibly deprecated"""
    return False


def test_file(value):
    """Possibly deprecated"""
    return False



def test_int(value):
    """Check if an input value is an integer, not sure if this is used, possibly deprecated"""
    try:
        floatval = float(value)
        intval = int(value)
        if intval == floatval and "." not in unicode(value):
            return True
    except  ValueError:
        return False

def test_number(value):
    """Check if an input value is an number, not sure if this is used, possibly deprecated"""
    try:
        floatval = float(value)
        return True
    except  ValueError:
        return False


def test_stringdate(value):
    """Check if an input value is an date, not sure if this is used, possibly deprecated"""
    try:
        curated_value = dateutil.parser.parse(unicode(value)).strftime("%Y-%m-%d")
        return curated_value
    except ValueError:
        return False

def test_percentage(value):
    """Check if an input value is an percentage, not sure if this is used, possibly deprecated"""
    result = test_number(value)
    if result:
        if float(val) > 0 and float(val) < 100:
            return True
    return False




class PinnedCustomField(TimeStampedModel):
    """PinnedCustomField is the model 
    which stores all of the information about 
    a particular field in a project"""
    TEXT = "text"
    TEXTAREA = "textarea"
    UISELECT = "uiselect"
    INTEGER = "integer"
    NUMBER = "number"
    RADIOS = "radios"
    UISELECTTAG = "uiselecttag"
    UISELECTTAGS = "uiselecttags"
    CHECKBOXES = "checkboxes"
    PERCENTAGE = "percentage"
    DECIMAL = "decimal"
    BOOLEAN = "boolean"
    DATE = "date"
    IMAGE = "imghref"
    LINK = "href"
    FILE_ATTACHMENT = "file"


    def pandas_converter(self, field_length, pandas_dtype):
        """Convert a variable from its pandas datatype, not sure if used podssibly deprecated"""
        dtype = str(pandas_dtype)
        if dtype == "int64":
            return self.INTEGER
        if dtype == "float64":
            return self.NUMBER
        if dtype == "object":
            if field_length > 100:
                return self.TEXTAREA
            else:
                return self.TEXT
        return self.TEXT

    FIELD_TYPE_CHOICES = OrderedDict((
        (TEXT, {"name": "Short text field", "display_form":{"type" :"copyfield"} ,"data": {"className": "htCenter htMiddle ",
         "type": "string", "icon": "<span class ='glyphicon glyphicon-font'></span>", "renderer_named" : "defaultCustomFieldRenderer" }, "test_datatype" : test_string}
         ),
        ("char", {
            "name": "Short text field", "display_form":{"type" :"copyfield"}, "data": {"className": "htCenter htMiddle ","type": "string", "renderer_named" : "defaultCustomFieldRenderer"}, "test_datatype" : test_string}),
        # (TEXTAREA, {"name": "Full text", "data": {
        #  "icon": "<span class ='glyphicon glyphicon-font'></span>", "type": "string", "format": "textarea"}, "test_datatype" : test_string}),
        (TEXTAREA, {"name": "Full text", "display_form":{"type" :"copyfield"}, "data": { "className": "htCenter htMiddle ","renderer_named" : "defaultCustomFieldRenderer",
         "icon": "<span class ='glyphicon glyphicon-font'></span>", "type": "string", "format": "ckeditor", "ckeditor":{ 'config.toolbar_Standard': [ { 'name': 'basicstyles', 'items' : [ 'Italic','Underline','Strike','Subscript','Superscript','-','RemoveFormat' ] }, ] },}, "test_datatype" : test_string}),
        (UISELECT, {"name": "Choice field", "display_form":{"type" :"copyfield"}, 
            "data": {"className": "htCenter htMiddle ",
         "type": "string", "format": "filtereddropdown", 
         "renderer_named" : "defaultCustomFieldRenderer",
         "options":{
                          "tagging" : False,
                          "fetchDataEventName" : "openedTaggingDropdown",
                          "dataArrivesEventName" : "autoCompleteTaggingData",
                          "multiple" : False,
                          "staticItems" : [] 
                },

         }, "test_datatype" : test_string}),
        (INTEGER, {"name": "Integer field", "display_form":{"type" :"copyfield"}, "data": {"className": "htCenter htMiddle ",
         "icon": "<span class ='glyphicon glyphicon-stats'></span>", "type": "integer", "renderer_named" : "defaultCustomFieldRenderer"}, "test_datatype" : test_int}),
        (NUMBER, {"name": "Decimal field", "display_form":{"type" :"copyfield"}, "data": {"className": "htCenter htMiddle ",
         "icon": "<span class ='glyphicon glyphicon-sound-5-1'></span>", "type": "number", "renderer_named" : "defaultCustomFieldRenderer"}, "test_datatype": test_number}),
        (UISELECTTAG, {"name": "Choice allowing create", "display_form":{"type" :"copyfield"},
         "data":  {"className": "htCenter htMiddle ",
         "icon": "<span class ='glyphicon glyphicon-tag'></span>",
          "type": "string", 
          "format": "filtereddropdown", 
          "options":{
                          "tagging" : True,
                          "fetchDataEventName" : "openedTaggingDropdown",
                          "dataArrivesEventName" : "autoCompleteTaggingData",
                          "multiple" : False,
                          "staticItems" : [] 
                }, "renderer_named" : "defaultCustomFieldRenderer"}, "test_datatype" : test_string}),
        (UISELECTTAGS, {"name": "Tags field allowing create", 
            "display_form":{"type" :"copyfield"}, 
            "data": 
            {"className": "htCenter htMiddle ",
            "renderer_named" : "defaultCustomFieldRenderer" 
            , "icon": "<span class ='glyphicon glyphicon-tags'></span>", 
            "type": "array", 
            "format": "filtereddropdown", 
            "items": {
                    "type": "string"
                  },
             "options":{
                          "tagging" : True,
                          "fetchDataEventName" : "openedTaggingDropdown",
                          "dataArrivesEventName" : "autoCompleteTaggingData",
                          "multiple" : True,
                          "staticItems" : [] 
                }
            }, "test_datatype" : test_string
            }),
        (PERCENTAGE, {"name": "Percentage field", "display_form":{"type" :"copyfield"}, "data": {"className": "htCenter htMiddle ","renderer_named" : "defaultCustomFieldRenderer",
         "icon": "<span class ='glyphicon'>%</span>", "type": "number", "maximum": 100.0, "minimum": 0.1}, "test_datatype": test_percentage}),
        (RADIOS, {"name": "Radio buttons", "display_form":{"type" :"copyfield"}, "data": {"className": "htCenter htMiddle ","renderer_named" : "defaultCustomFieldRenderer",
         "icon": "<span class ='glyphicon'>%</span>", "type": "string", "format":"radios"}, "test_datatype": test_string}),
        (DATE,  {"name": "Date Field", "display_form":{"type" :"copyfield"}, "data": {"className": "htCenter htMiddle ", "renderer_named" : "defaultCustomFieldRenderer",
         "icon": "<span class ='glyphicon glyphicon-calendar'></span>", "type": "string",   "format": "date"}, "test_datatype": test_stringdate}),
        (LINK, {"name": "Link to server or external", "display_form":{"type" :"copyfield"}, "data": {"className": "htCenter htMiddle ","renderer_named" : "defaultCustomFieldRenderer", "format": "href", "type":
                                                               "string", "icon": "<span class ='glyphicon glyphicon glyphicon-new-window'></span>"}, "test_datatype" : test_string}),
        (IMAGE, {"name": "Image link to embed", "display_form":{"type" :"copyfield"}, "data": {"className": "htCenter htMiddle ","renderer_named" : "defaultCustomFieldRenderer" , "format": "imghref", "type":
                                                         "string", "icon": "<span class ='glyphicon glyphicon glyphicon-picture'></span>"},  "test_datatype" : test_string}),
        (DECIMAL, {"name": "Decimal field", "display_form":{"type" :"copyfield"}, "data": {"className": "htCenter htMiddle ",
         "renderer_named" : "defaultCustomFieldRenderer", "icon": "<span class ='glyphicon'>3.1</span>", "type": "number"},  "test_datatype" : test_number}),
        (BOOLEAN, {"name": "checkbox", "display_form":{"type" :"copyfield"}, "data": {"className": "htCenter htMiddle ","renderer_named" : "defaultCustomFieldRenderer",
         "icon": "<span class ='glyphicon'>3.1</span>", "type": "boolean"},"test_datatype" : test_bool }),
        ("related", {"name": "TEST", "display_form":{"type" :"copyfield"}, "data": {"className": "htCenter htMiddle ", "renderer_named" : "defaultCustomFieldRenderer",
         "icon": "<span class ='glyphicon'>3.1</span>", "type": "string"}, "test_datatype" : test_string}),
        (FILE_ATTACHMENT, {"name": "File Upload", "display_form":{"type" :"attachmentlist"}, "data": { "className": "htCenter htMiddle ",
         "icon": "<span class ='glyphicon glyphicon-paperclip'></span>", "type": "object", "format": "file_upload", "renderer_named" : "fileUploadRenderer"}, "test_datatype" : test_file})

    ))

    field_key = models.CharField(max_length=500,  default="", help_text="field key value, not currently used perhaps deprecated")
    name = models.CharField(max_length=500, null=False, blank=False, help_text="Name of the field in the project, may contain any character apart from slashes and dots")
    description = models.CharField(
        max_length=1024, blank=True, null=True, default="", help_text="Description of the field to be displayed in the angular schema form when editing data")
    custom_field_config = models.ForeignKey(
        "cbh_core_model.CustomFieldConfig", related_name='pinned_custom_field', default=None, null=True, blank=True, help_text="the custom field config object this field is part of")
    required = models.BooleanField(default=False, help_text="Whetehr the field is required")
    field_type = models.CharField(default="char", choices=(
        (name, value["name"]) for name, value in FIELD_TYPE_CHOICES.items()), max_length=15, help_text="The data type of the data to be added in the front end forms for this field" )
    allowed_values = models.CharField(
        max_length=1024, blank=True, null=True, default="", help_text="What values are allowed for this field")
    position = models.PositiveSmallIntegerField(help_text="Auto-filled field that says what order the fields in a custom field config should be displayed int he form")
    default = models.CharField(max_length=500, default="", blank=True, help_text="The default value of this field to be applied when adding data via the angular schema form")

    open_or_restricted = models.CharField(max_length=20, default=OPEN, choices=RESTRICTION_CHOICES,  help_text="Whether to open up this field to people who only have viewer rights on the project")

    def validate_field(self, value):
        """Data type testing for fields in the custom fields of a compound batch (possibly deprecated or unfinished"""
        if not value and self.required:
            return False
        else:
            func = self.FIELD_TYPE_CHOICES[self.field_type]["test_datatype"]
            return func(value)


    @cached_property
    def get_items_simple(self):
        """List the allowed values for a particular field"""
        items = [item.strip()
                 for item in self.allowed_values.split(",") if item.strip()]
        setitems = sorted(list(set(items)))
        testdata = [{"doc_count": 0, "key": item.strip()}
                    for item in setitems if item]
        return testdata

    @cached_property
    def get_space_replaced_name(self):
        """Return the name of the field in a space replaced way for ealsticsearch (deprecated)"""
        return self.name.replace(u" ", u"__space__")

    def __unicode__(self):
        """Unicode representation of the field"""
        return "%s  %s  %s" % (self.name, self.field_key,  self.field_type)

    @cached_property
    def field_values(obj):
        """Pull out the data that is required when compiling the angular schema form JSON for a particular field"""
        data = copy(obj.FIELD_TYPE_CHOICES[obj.field_type]["data"])

        data["title"] = obj.name
        data["placeholder"] = obj.description
        
        display_form = copy(obj.FIELD_TYPE_CHOICES[obj.field_type]["display_form"])
        form = {}
        form["knownBy"] = obj.name
        form["data"] = "custom_fields.%s" % obj.name
        form["position"] = obj.position
        form["key"] = obj.name
        display_form["key"] = obj.get_space_replaced_name
        form["title"] = obj.name

        form["description"] = obj.description
        form["disableSuccessState"] = True
        form["feedback"] = True
        searchitems = []
        if obj.default:
            data['default'] = obj.default
        else:
            data['default'] = ""

        if data["type"] == "array":
            if obj.default:
                data['default'] = obj.default.split(",")
            else:
                data['default'] = []

            if obj.required:
                form["minLength"] = 1

        if "filtereddropdown" in data.get("format", ""):
            form["type"] = "filtereddropdown"
            form["placeholder"] = "Choose..."
            data['options']['staticItems'] = obj.get_items_simple

        if "radios" in data.get("format", ""):
            form["type"] = "radios"
            data['enum'] = [value["key"] for value in obj.get_items_simple]
            form['titleMap'] = [{"name": value["key"], "value": value["key"]} for value in obj.get_items_simple]


        if data.get("format", False) == "file_upload":
            #will need to alter the init method here (so it's no longer using dataoverviewctrl)
            #also the success method to add the projectid? Not sure what that means really
            form["uploadOptions"] = { "modal": {
                                          'title': 'Modal Title specified here',
                                          'flow': {
                                            'dropEnabled': False,
                                            'imageOnly': False,
                                            'init': {},
                                            'success': 'success(file, formkey)',
                                            'removeFile': 'removeFile(formkey, index, url)',
                                            'imageFunction': 'fetchImage(url)',
                                            'sizeCheck': 'sizeCheck(file, formkey)'
                                          }
                                        }
                                      }
            data['default'] = {"attachments" : []}
            form["default"] = {"attachments" : []}



        #add config options for ckeditor
        if data.get("format", False) == "ckeditor":
            form['ckeditor'] = { 'toolbar': [
                                                { 'name': 'clipboard', 'items': [ 'Cut', 'Copy', 'Paste', 'PasteText', 'PasteFromWord', '-', 'Undo', 'Redo' ] },
                                                { 'name': 'editing', 'items': [ 'Scayt' ] },
                                                { 'name': 'links', 'items': [ 'Link', 'Unlink', 'Anchor' ] },
                                                { 'name': 'tools', 'items': [ 'Maximize' ] },
                                                '/',
                                                { 'name': 'basicstyles', 'items' : [ 'Bold','Italic','Underline','Strike','Subscript','Superscript','-','RemoveFormat' ] },
                                                { 'name': 'insert', 'items' : [ 'HorizontalRule','SpecialChar','PageBreak' ] },
                                                { 'name': 'paragraph', 'items' : [ 'NumberedList','BulletedList','-','Outdent','Indent','-','JustifyLeft','JustifyCenter','JustifyRight','JustifyBlock'] },
                                                { 'name': 'styles', 'items' : [ 'Styles', 'Format' ] },
                                                { 'name': 'about', 'items' : [ 'About' ] }
                                            ] }

        if data.get("format", False) == obj.DATE:
            maxdate = time.strftime("%Y-%m-%d")
            form.update({
                "minDate": "2000-01-01",
                "maxDate": maxdate,
                'type': 'datepicker',
                "format": "yyyy-mm-dd",
                'pickadate': {
                    'selectYears': True,
                    'selectMonths': True,
                },
            })

        else:
            for item in ["options"]:
                stuff = data.pop(item, None)
                if stuff:
                    form[item] = stuff
        return (deepcopy(data) , deepcopy(form), deepcopy(display_form))

    class Meta:
        """Ordering by default in the model for use with tastypie related resource which does not give a way to order the related field values"""
        ordering = ['position']
        get_latest_by = 'created'


class Invitation(TimeStampedModel):
    """Invitation model which saves the fact that an invitation has been sent to a given user"""
    email = models.CharField(max_length=100)
    created_by = models.ForeignKey("auth.User")
    first_name = models.TextField(default="", null=True, blank=True, )
    last_name = models.TextField(default="", null=True, blank=True, )
    projects = models.ManyToManyField("cbh_core_model.Project", blank=True)

    def __unicode__(self):
        return self.name




class CBHFlowFile(models.Model):
    """
    A file upload through Flow.js
    """
    STATE_UPLOADING = 1
    STATE_COMPLETED = 2
    STATE_UPLOAD_ERROR = 3

    STATE_CHOICES = [
        (STATE_UPLOADING, "Uploading"),
        (STATE_COMPLETED, "Completed"),
        (STATE_UPLOAD_ERROR, "Upload Error")
    ]

    # identification and file details
    identifier = models.SlugField(max_length=255, unique=True, db_index=True, help_text="String identifier built using the name and session id of the file for security so other sessions cannot access the file")
    original_filename = models.CharField(max_length=200, help_text="The original filename of the file")
    total_size = models.IntegerField(default=0, help_text="Size in bytes")
    total_chunks = models.IntegerField(default=0, help_text="Number of chunks the file was split up into to upload")

    # current state
    total_chunks_uploaded = models.IntegerField(default=0, help_text="Number of chunks uploaded so far")
    state = models.IntegerField(choices=STATE_CHOICES, default=STATE_UPLOADING, help_text="Current status of the upload")
    created = models.DateTimeField(auto_now_add=True, help_text="Date the upload was created")
    updated = models.DateField(auto_now=True, help_text="Date the upload was updated")
    project = models.ForeignKey("cbh_core_model.Project", help_text="Project that the uploaded file is associated with")

    def __unicode__(self):
        """Unicode representation of the file"""
        return self.identifier

    def update(self):
        """Check the status of the uploading chunks"""
        self.total_chunks_uploaded = self.chunks.count()
        super(CBHFlowFile, self).save()
        self.join_chunks()

    @property
    def extension(self):
        """
        Return the extension of the upload file
        """
        name, ext = os.path.splitext(self.original_filename)
        return ext

    @property
    def filename(self):
        """
        Return the unique filename generated based on the identifier
        """
        return '%s%s' % (self.identifier, self.extension)

    @property
    def file(self):
        """
        Return the uploaded file
        """
        if self.state == self.STATE_COMPLETED:
            return default_storage.open(self.path)
        return None

    @property
    def full_path(self):
        """
        Return the full path of the file uploaded
        """
        return os.path.join(settings.MEDIA_ROOT, self.path)

    @property
    def path(self):
        """
        Return the path of the file uploaded
        """
        return os.path.join(FLOWJS_PATH, self.filename)

    def get_chunk_filename(self, number):
        """
        Return the filename of the chunk based on the identifier and chunk number
        """
        return '%s-%s.tmp' % (self.identifier, number)

    def join_chunks(self):
        """
        Join all the chunks in one file
        """
        if self.state == self.STATE_UPLOADING and self.total_chunks_uploaded == self.total_chunks:

            # create file and write chunks in the right order
            temp_file = open(self.full_path, "wb")
            for chunk in self.chunks.all():
                chunk_bytes = chunk.file.read()
                temp_file.write(chunk_bytes)
            temp_file.close()

            # set state as completed
            self.state = self.STATE_COMPLETED
            super(CBHFlowFile, self).save()

            # delete chunks automatically if is activated in settings
            if FLOWJS_AUTO_DELETE_CHUNKS:
                self.chunks.all().delete()

    def is_valid_session(self, session):
        """
        Check if a session id is the same that uploaded the file
        """
        return self.identifier.startswith(session)


class CBHFlowFileChunk(models.Model):
    """
    A chunk is part of the file uploaded
    """
    class Meta:
        ordering = ['number']

    # identification and file details
    parent = models.ForeignKey(CBHFlowFile, related_name="chunks", )
    file = models.FileField(max_length=255, upload_to=chunk_upload_to)
    number = models.IntegerField()
    created_at = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return self.filename

    @property
    def filename(self):
        return self.parent.get_chunk_filename(self.number)

    def save(self, *args, **kwargs):
        super(CBHFlowFileChunk, self).save(*args, **kwargs)
        self.parent.update()


@receiver(pre_delete, sender=CBHFlowFile)
def flow_file_delete(sender, instance, **kwargs):
    """
    Remove files on delete if is activated in settings
    """
    if FLOWJS_REMOVE_FILES_ON_DELETE:
        try:
            default_storage.delete(instance.path)
        except NotImplementedError:
            pass


@receiver(pre_delete, sender=CBHFlowFileChunk)
def flow_file_chunk_delete(sender, instance, **kwargs):
    """
    Remove file when chunk is deleted
    """
    instance.file.delete(False)



def print_name(sender, instance, **kwargs):
    instance.project_key = slugify(instance.name)

    try:
        print sender.objects.get(sender.id).__dict__
    except:
        print "new"

pre_save.connect(print_name, sender=Project, dispatch_uid="proj_key1")
pre_save.connect(print_name, sender=CustomFieldConfig, dispatch_uid="proj_key2")
pre_save.connect(print_name, sender=PinnedCustomField, dispatch_uid="proj_key3")

