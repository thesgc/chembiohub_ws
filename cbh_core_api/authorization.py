"""Main authorization module for ChemBio Hub platform API"""

from tastypie.authorization import Authorization
from tastypie.exceptions import Unauthorized
import logging
logger = logging.getLogger(__name__)
logger_debug = logging.getLogger(__name__)
from cbh_core_model.models import Project, get_all_project_ids_for_user_perms, get_all_project_ids_for_user, RESTRICTED, get_projects_where_fields_restricted, PERMISSION_CODENAME_SEPARATOR

def viewer_projects(user):
    """Return the project ids which a particular user has viewer rights on"""
    pids = get_all_project_ids_for_user(user, ["viewer","editor", "owner"])
    return pids

def editor_projects(user):
    """Return the project ids which a particular user has editor rights on"""
    pids = get_all_project_ids_for_user(user, ["editor", "owner"])
    return pids

def owner_projects(user):
    """Return the project ids which a particular user has owner rights on"""
    pids = get_all_project_ids_for_user(user, ["owner"])
    return pids



class ProjectPermissionAuthorization(Authorization):
    """Authorization for the update of the permissions against a particular project"""
    def login_checks(self, request):
        """standard login checks to see if user is logged in"""
        if not hasattr(request, 'user'):
            print "no_logged_in_user"
            raise Unauthorized("no_logged_in_user")
        if not request.user.is_authenticated():
            raise Unauthorized("no_logged_in")

    def create_detail(self, object_list, bundle):
        """We do not yet support creating of new permissions by API so return false"""
        return False


    def update_list(self, object_list, bundle, for_list=True):
        """
        If any particular item cannot be updated then raise an exception
        """
        self.login_checks(bundle.request)
        for permissionbundle in bundle.data["objects"]:
            self.update_detail(object_list, permissionbundle)
        return object_list

    def update_detail(self, object_list, bundle, for_list=False):
        """
        Function is called by the patch detail method of the project permission resource
        Ensure that the person updating a project has owner rights
        Ensure that the owner is not removing themselves
        Ensure that when a particular user is 
        being given ownership of a project that they 
        have the add project permission
        """
        if not for_list:
            self.login_checks(bundle.request)
        try:
            entity_id, perm_code  = tuple(bundle.obj.codename.split(PERMISSION_CODENAME_SEPARATOR))
            #check that the user changing the permissions for a project has owner rights on that project
            

            if entity_id.isdigit():
                if int(entity_id) not in owner_projects(bundle.request.user):
                    raise Unauthorized("You are not an owner of the project with id %s so you cannot change its permissions" % entity_id)
                #If we are updating the owner of a certain project, check that every person in the list being added as an owner has the "add project" privilege
                if perm_code == "owner":
                    owner_user_is_in_new_list = False
                    for userbun in bundle.data["users"]:
                        if userbun.obj.id == bundle.request.user.id:
                            owner_user_is_in_new_list = True
                        if not userbun.obj.has_perm("cbh_core_model.add_project"):
                            raise Unauthorized("You cannot add owner permissions to this user")
                    if bundle.request.user.is_superuser is False:
                        if owner_user_is_in_new_list is False:
                            raise Unauthorized("You may not remove yourself as an owner of a project")

        except ValueError:
            raise Unauthorized("You are not permitted to update the permission with codename %s via the API")
        return True



class InviteAuthorization(Authorization):
    """Authorization of invitations to ChemBio Hub Platform"""
    def login_checks(self, request, model_klass, perms=None):
        """standard login checks to see if user is logged in"""

        if not hasattr(request, 'user'):
            print "no_logged_in_user"
            raise Unauthorized("no_logged_in_user")
        if not request.user.is_authenticated():
            raise Unauthorized("no_logged_in")


    def create_list(self, object_list, bundle):
        return []


    def create_detail(self, object_list, bundle):
        """Checks that the person doing the invting is an editor
        This function is called during the obj_create function in tastypie
        """
        self.login_checks(bundle.request, bundle.obj.__class__)
        pids = editor_projects(bundle.request.user)
        for project in bundle.data["projects_selected"]:
            if(project["id"] not in pids):
                raise Unauthorized("Not authorized to invite to this project, you must have editor status")
        return True


    def update_detail(self, object_list, bundle):
        """There is not yet a use case for updating a 
        single invite record"""
        raise Unauthorized("not authroized for to update")







class ProjectListAuthorization(Authorization):

    """
    Permissions implementation for the Projects API
    """

    

    def login_checks(self, request, model_klass):
        """standard login checks to see if user is logged in"""

        if not hasattr(request, 'user'):
            print "no_logged_in_user"
            raise Unauthorized("no_logged_in_user")


    def list_checks(self, request, model_klass, object_list):
        """Filter down a queryset of projects to give only the ones the user has permission to see"""
        perms = request.user.get_all_permissions()
        pids = viewer_projects(request.user)
        self.login_checks(request,  model_klass, )

        return object_list.filter(pk__in=pids)


    def alter_project_data_for_permissions(self, bundle, request):
        """Add appropriate flags and tidy up the restricted fields in accordance with
        the permissions that a particular user has"""
        edit_projects = editor_projects(request.user)
        own_projects = owner_projects(request.user)
        restricted_and_unrestricted_projects = get_projects_where_fields_restricted(request.user)

        if isinstance(bundle, dict):
            for bun in bundle['objects']:
                bun.data['editor'] = bun.obj.id in edit_projects
                bun.data["owner"] = bun.obj.id in own_projects
                self.alter_bundle_for_user_custom_field_restrictions(bun, restricted_and_unrestricted_projects)
        else:
            bundle.data['editor'] = bundle.obj.id in edit_projects
            bundle.data["owner"] = bundle.obj.id in own_projects






    def alter_bundle_for_user_custom_field_restrictions(self, bundle, restricted_and_unrestricted_projects):
        """Post serialization modification to the list of fields based on the field permissions
        We add flags to the project objects to be used in the front end to see if they are an editor or viewer 
        and restrict functionality accordingly (naturally back end validation is also done)
        We remove all of the restriucted fields so that the user never even knows they exist.
        Given that the user will be a viewer of the project not an editor, thisb does not affect editing of data
        """

        if bundle.data["id"] in restricted_and_unrestricted_projects[RESTRICTED]:
            new_fields = []
            for field in bundle.data["custom_field_config"].data["project_data_fields"]:
                if field.data["open_or_restricted"] == RESTRICTED:
                    #This is a restricted field and the user's access is restricted therefore block them
                    bundle.data["users_restricted_fields"].append(field["name"])
                else:
                    new_fields.append(field)
            bundle.data["custom_field_config"].data["project_data_fields"] = new_fields
            


    def read_list(self, object_list, bundle):
        """Filter down the 'get_list' objects that are serialized based on the user's permissions"""
        return self.list_checks(bundle.request, bundle.obj.__class__, object_list)

    def create_list(self, object_list, bundle):
        '''Creating lists is not allowed'''
        raise Unauthorized("Creating lists of projects not supported")
        
    def create_detail(self, object_list, bundle):
        '''Creating projects is allowed for all logged in users'''
        self.login_checks(bundle.request, bundle.obj.__class__)
        return bundle.request.user.has_perm("cbh_core_model.add_project")
        

    def update_list(self, object_list, bundle):
        '''Only owners of projects are allowed to update them'''
        return []

    def update_detail(self, object_list, bundle):
        '''Only owners of projects are allowed to update them'''
        self.login_checks(bundle.request, bundle.obj.__class__)
        pids = owner_projects(bundle.request.user)
        if bundle.obj.id in pids:
            return True
        raise Unauthorized("Not authorized to update project")




class ProjectAuthorization(Authorization):

    """
    A Permission class to be used against any object which is linked 
    by foreign key relationship to the project object
    """

    def login_checks(self, request, model_klass, perms=None):
        """standard login checks to see if user is logged in"""

        if not hasattr(request, 'user'):
            print "no_logged_in_user"
            raise Unauthorized("no_logged_in_user")
        if not request.user.is_authenticated():
            raise Unauthorized("no_logged_in")

 

    def project_ids(self, request ):
        """Commonly used function to get back the projects a user can view"""
        self.login_checks( request, None)
        pids = viewer_projects(request.user)
        return pids

    def create_list(self, object_list, bundle):
        """Not implemented"""
        return []

    def read_detail(self, object_list, bundle):
        """probably deprecated"""
        self.login_checks(bundle.request, bundle.obj.__class__)
        pids = viewer_projects(bundle.request.user)
        if bundle.obj.project.id in pids:
            return True
        else:
            raise Unauthorized("not authroized for project")

    def update_list(self, object_list, bundle):
        """Not used"""

        return []

    def create_detail(self, object_list, bundle):
        """Check for a pk in the project object in order to see if the user is allowed to create a compoundbatch or other object in theat project"""
        self.login_checks(bundle.request, bundle.obj.__class__)
        pids = editor_projects(bundle.request.user)
        id = None
        if hasattr(bundle.data["project"], "id"):
            if bundle.data["project"].id in pids:
                return True
        elif hasattr(bundle.data["project"], "get"):
            if bundle.data["project"].get("pk", 0) in pids:
                return True

        raise Unauthorized("not authroized for project")
        # return self.base_checks(bundle.request, bundle.obj.__class__,
        # bundle.data, ["editor",])

    def update_detail(self, object_list, bundle):
        """Check if a user is allowed to update a specific record"""
        self.login_checks(bundle.request, bundle.obj.__class__)
        pids = editor_projects(bundle.request.user)
        if bundle.obj.project.id in pids:
            return True

        raise Unauthorized("not authroized for project")

    def read_list(self, object_list, bundle):
        """Read list is now implemented in elasticsearch instead"""
        return object_list

    def check_if_field_restricted(self, field_path, project_ids_requested, tabular_data_schema):
        """Here we take the field path value which points to a custom field and look it up in the global schema to projects_selected
        if it is restricted for any of the projects in the list that the user has chosen for this request"""

        #because this schema has been requested on a per-user basis, the fields which are restricted have already been taken into account and removed
        #Therefore the project specific schema contains projects for which this field name is not restricted

        proj_specific_schema = tabular_data_schema["schema"][field_path]

        return [pidr  for pidr in project_ids_requested if pidr in proj_specific_schema]


    def remove_restricted_fields_from_custom_fields_to_render(self, batch_dicts, user_restricted_fieldnames):
        """Given a list of the field names to be removed on a per project basis - remove them based on the project id in the batch dictionary"""

        for fieldname in user_restricted_fieldnames[batch_dict["projectfull"]["id"]]:
                if fieldname in batch_dict["custom_fields"]:
                    del batch_dict["custom_fields"][fieldname]
        return batch_dict
