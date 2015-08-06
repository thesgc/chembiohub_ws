from django.db import models, connection

from solo.models import SingletonModel


def get_all_hstore_values(table,column, key, is_list=False, extra_where=" True"):
    '''Using an hstore query from the reference here
    http://www.youlikeprogramming.com/2012/06/mastering-the-postgresql-hstore-date-type/
    where project_id = {{project_id}}
    '''
    cursor = connection.cursor()
    sql = "SELECT DISTINCT {column} -> '{key}' FROM {table} where {column} -> '{key}' != '' and {extra_where};".format( **{"table":table, "key":key, "column": column,  "extra_where": extra_where})
    cursor.execute(sql)
    mytuple = cursor.fetchall()
    items = []
    data = [d[0] for d in mytuple]
    for d in data:
        if is_list:
            try:
                for elem in json.loads(d):
                    items.append(elem)
            except ValueError:
                items.append(d)
            except TypeError:
                items.append(d)
        else:
            items.append(d)
    return items

PROJECT_PERMISSIONS = (("viewer","Can View"),
                        ("editor","Can edit or add batches"),
                        ( "admin", "Can assign permissions"))

class ProjectPermissionManager(models.Manager):

    def sync_all_permissions(self):
        for perm in self.all():
            perm.sync_permissions()

    def get_user_permission(self,project_id, user, codenames, perms=None):
        '''Check the given users' permissions against a list of codenames for a project id'''
        if not perms:
            perms = user.get_all_permissions()
        codes = ["%d.%s" % (project_id, codename) for codename in codenames]
        matched = list(perms.intersection(codes))
        if len(matched) > 0:
            return True
        return False

    



class ProjectPermissionMixin(models.Model):
    '''The aim of this mixin is to create a permission content type and a permission model for a given project
    It allows for pruning the contnet types once the model is changed
    '''
    

    objects = ProjectPermissionManager()
    
    def get_project_key(self):
        return str(self.pk)

    def sync_permissions(self):
        '''first we delete the existing permissions that are not labelled in the model'''
        ct , created = ContentType.objects.get_or_create(app_label=self.get_project_key(), model=self, name=self.name)
        deleteable_permissions = Permission.objects.filter(content_type_id=ct.pk).exclude(codename__in=[perm[0] for perm in PROJECT_PERMISSIONS])
        deleteable_permissions.delete()
        for perm in PROJECT_PERMISSIONS:
            pm = Permission.objects.get_or_create(content_type_id=ct.id,codename=perm[0],name=perm[1])



    def get_contenttype_for_instance(self):
        ct = ContentType.objects.get(app_label=self.get_project_key(), model=self)
        return ct

    def delete_all_instance_permissions(self):
        '''for the pre delete signal'''
        deleteable_permissions = Permission.objects.filter(content_type_id=self.get_contenttype_for_instance().id)
        deleteable_permissions.delete()
                

    def get_instance_permission_by_codename(self, codename):
        pm = Permission.objects.get(codename=codename, content_type_id=self.get_contenttype_for_instance().id)
        return pm


    def _add_instance_permissions_to_user_or_group(self, group_or_user, codename):
        if type(group_or_user) == Group:
            group_or_user.permissions.add(self.get_instance_permission_by_codename(codename))
        if type(group_or_user) == User:
            group_or_user.user_permissions.add(self.get_instance_permission_by_codename(codename))

    def make_editor(self,group_or_user):
        self._add_instance_permissions_to_user_or_group(group_or_user, "editor")


    def make_viewer(self,group_or_user):
        self._add_instance_permissions_to_user_or_group(group_or_user, "viewer")


    def make_admin(self,group_or_user):
        self._add_instance_permissions_to_user_or_group(group_or_user, "admin")



    class Meta:
       
        abstract = True
        



class ProjectType(TimeStampedModel):
    ''' Allows configuration of parts of the app on a per project basis - initially will be used to separate out compound and inventory projects '''
    name = models.CharField(max_length=100, db_index=True, null=True, blank=True, default=None)
    show_compounds = models.BooleanField(default=True)

    def __unicode__(self):
        return self.name






class Project(TimeStampedModel, ProjectPermissionMixin):
    ''' Project is a holder for moleculedictionary objects and for batches'''
    name = models.CharField(max_length=50, db_index=True, null=True, blank=True, default=None)
    project_key = models.SlugField(max_length=50, db_index=True, null=True, blank=True, default=None, unique=True)
    created_by = models.ForeignKey("auth.User")
    custom_field_config = models.ForeignKey("cbh_chembl_model_extension.CustomFieldConfig", related_name="project",null=True, blank=True, default=None, )
    project_type = models.ForeignKey(ProjectType,null=True, blank=True, default=None)
    is_default = models.BooleanField(default=False)
    #allow configuration of project types 
    #project_type = models.
    
    class Meta:
        get_latest_by = 'created'

    def __unicode__(self):
        return self.name

    @models.permalink
    def get_absolute_url(self):
        return {'post_slug': self.project_key}



def sync_permissions(sender, instance, created, **kwargs):
    '''After saving the project make sure it has entries in the permissions table'''
    if created is True:
        instance.sync_permissions()
        instance.make_editor(instance.created_by)

post_save.connect(sync_permissions, sender=Project, dispatch_uid="proj_perms")





class CustomFieldConfig(TimeStampedModel):
    name = models.CharField(unique=True, max_length=50)
    created_by = models.ForeignKey("auth.User")
    schemaform = models.TextField(default = "", null=True, blank=True, )
    def __unicode__(self):
        return self.name

    

class CBHCompoundMultipleBatch(TimeStampedModel):
    '''Holds a list of batches'''
    created_by = models.CharField(max_length=50, db_index=True, null=True, blank=True, default=None)
    project = models.ForeignKey(Project, null=True, blank=True, default=None)
    uploaded_data = PickledObjectField()
    uploaded_file = models.OneToOneField(FlowFile, null=True, blank=True, default=None)
    saved = models.BooleanField(default=False)
    #batches = models.ForeignKey(CBHCompoundBatch, null=True, default=None)

class SkinningConfig(SingletonModel):
    '''Holds information about custom system messages and other customisable elements'''
    #created_by = models.ForeignKey("auth.User")
    instance_alias = models.CharField(max_length=50, null=True, blank=False, default='ChemReg')
    project_alias = models.CharField(max_length=50, null=True, blank=False, default='project')
    result_alias = models.CharField(max_length=50,null=True, blank=False, default='result')

    def __unicode__(self):
        return u"Skinning Configuration"

    class Meta:
        verbose_name = "Skinning Configuration"
    #we can eventually use this to specify different chem sketching tools
    









class PinnedCustomField(TimeStampedModel):
    TEXT = "text"
    TEXTAREA = "textarea"
    UISELECT = "uiselect"
    INTEGER  = "integer"
    NUMBER = "number"
    UISELECTTAG  = "uiselecttag"
    UISELECTTAGS = "uiselecttags"
    CHECKBOXES = "checkboxes"
    PERCENTAGE = "percentage"
    DATE = "date"
    IMAGE = "imghref"
    LINK = "href"


    FIELD_TYPE_CHOICES = {
                            "char" : {"name" : "Short text field", "data": { "type": "string"}},

                            TEXT : {"name" : "Short text field", "data": { "type": "string" ,"icon":"<span class ='glyphicon glyphicon-font'></span>" }},
                            TEXTAREA: {"name" :"Full text", "data": { "icon":"<span class ='glyphicon glyphicon-font'></span>","type": "string" , "format" : "textarea"}},
                            UISELECT: {"name" :"Choice field", "data": { "type": "string" , "format" : "uiselect"}},
                            INTEGER: {"name" :"Integer field", "data": { "icon":"<span class ='glyphicon glyphicon-stats'></span>" ,"type": "integer"}},
                            NUMBER: {"name" :"Decimal field", "data": { "icon":"<span class ='glyphicon glyphicon-sound-5-1'></span>","type": "number"}},
                            UISELECTTAG: {"name" : "Choice allowing create", "data":  { "icon":"<span class ='glyphicon glyphicon-tag'></span>", "type": "string", "format" : "uiselect"}},
                            UISELECTTAGS: {"name" : "Tags field allowing create" , "data": { "icon":"<span class ='glyphicon glyphicon-tags'></span>","type": "array", "format" : "uiselect", "options": {
                                      "tagging": "tagFunction" ,
                                      "taggingLabel": "(adding new)",
                                      "taggingTokens": "",
                                 }}},
                            PERCENTAGE: {"name" :"Percentage field", "data": { "icon":"<span class ='glyphicon'>%</span>", "type": "number", "maximum" : 100.0, "minimum": 0.1}},
                            DATE:  {"name": "Date Field" , "data":{"icon":"<span class ='glyphicon glyphicon-calendar'></span>","type": "string",   "format": "date"}},
                            LINK : {"name" : "Link to server or external", "data": { "format": "href", "type": "string" ,"icon":"<span class ='glyphicon glyphicon glyphicon-new-window'></span>" }},
                            IMAGE : {"name" : "Image link to embed", "data": {"format": "imghref", "type": "string" ,"icon":"<span class ='glyphicon glyphicon glyphicon-picture'></span>" }},
                        }


    field_key = models.CharField(max_length=50,  default="")
    name = models.CharField(max_length=50)
    description = models.CharField(max_length=1024, blank=True, null=True, default="")
    custom_field_config = models.ForeignKey("cbh_chembl_model_extension.CustomFieldConfig", related_name='pinned_custom_field')
    required = models.BooleanField(default=False)
    part_of_blinded_key = models.BooleanField(default=False, verbose_name="blind key")
    field_type = models.CharField(default="char", choices=((name, value["name"]) for name, value in FIELD_TYPE_CHOICES.items()), max_length=15, )
    allowed_values = models.CharField(max_length=1024, blank=True, null=True, default="")
    position = models.PositiveSmallIntegerField()

    def get_dropdown_list(self, projectKey):
        is_array = False
        if self.FIELD_TYPE_CHOICES[self.field_type]["data"]["type"] == "array":
            is_array=True
        db_items = get_all_hstore_values("cbh_chembl_model_extension_cbhcompoundbatch  inner join cbh_chembl_model_extension_project on cbh_chembl_model_extension_project.id = cbh_chembl_model_extension_cbhcompoundbatch.project_id ", 
            "custom_fields", 
            self.name, 
            is_list=is_array, 
            extra_where="cbh_chembl_model_extension_project.project_key ='%s'" % projectKey)
        return  [item for item in db_items]

    def get_allowed_items(self,projectKey):
        items = [item.strip() for item in self.allowed_values.split(",") if item.strip()]
        setitems = sorted(list(set(items + self.get_dropdown_list(projectKey))))
        testdata = [{"label" : item.strip(), "value": item.strip()} for item in setitems if item] 
        searchdata = [{"label" : "[%s] %s" % (self.name ,item.strip()), "value" : "%s|%s" % (self.name ,item.strip())} for item in setitems if item] 
        return (testdata, searchdata)

    class Meta:
        ordering = ['position']
        get_latest_by = 'created'


