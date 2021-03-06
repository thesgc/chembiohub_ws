from django.core.management.base import BaseCommand, CommandError
try: input = raw_input
except NameError: pass
import elasticsearch
import os
from django.contrib.auth import get_user_model
from django.test import RequestFactory
import shutil
from cbh_utils import elasticsearch_client
from cbh_core_api.tasks import remove_session_cached_projectlists

def backup_projects(projects_list, directory):
    from cbh_core_api.resources import ChemregProjectResource, ChemRegCustomFieldConfigResource
    cfcres = ChemRegCustomFieldConfigResource()
    projres = ChemregProjectResource()
    request_factory = RequestFactory()
    user = get_user_model().objects.filter(is_superuser=True)[0]
    for proj in projects_list:
        project_path = "%s/%d__%s/" % (directory, proj.id, proj.name )
        try:
            os.makedirs(project_path)
        except OSError:
            pass
        request = request_factory.get("/")
        request.user = user
        json_content = projres.get_detail(request, pk=proj.id).content
        path_to_write_json = project_path + "project.json"
        with open(path_to_write_json, "w") as text_file:
            text_file.write(json_content)

        request = request_factory.get("/")
        request.user = user
        request.GET = request.GET.copy()
        request.GET["format"] = "xlsx"
        resp = cfcres.get_detail(request, pk=proj.custom_field_config_id)
        path_to_write_xlsx = project_path + "fields.xlsx"
        with open(path_to_write_xlsx, 'wb') as f:
            f.write(resp._container[0])

        




def backup_compounds(projects_list, directory):
    from cbh_chem_api.resources import CBHCompoundBatchResource
    compoundres = CBHCompoundBatchResource()
    request_factory = RequestFactory()
    user = get_user_model().objects.filter(is_superuser=True)[0]
    
    for proj in projects_list:
        project_path = "%s/%d__%s/" % (directory, proj.id, proj.name )
        try:
            os.makedirs(project_path)
        except OSError:
            pass
        for format in ["sdf", "xlsx", "json"]:
            newrequest = request_factory.get("/")
            path_to_write = project_path + "compounds.%s" % format
            newrequest.user = user
            newrequest.GET = newrequest.GET.copy()
            newrequest.GET["format"] = format
            newrequest.GET["offset"] = 0
            newrequest.GET["limit"] = 100000
            newrequest.GET["pids"] = str(proj.id)
            file_resp = compoundres.get_list(newrequest)
            with open(path_to_write, 'wb') as f:
                f.write(file_resp._container[0])

def backup_filename(id, name):
    """get the backup filename for a file"""
    return str(id) + "__" + name


def backup_attachments(projects_list, directory):
    from cbh_core_model.models import CBHFlowFile

    for proj in projects_list:

        files_path = "%s/%d__%s/attachments/" % (directory, proj.id, proj.name )
        uploads_path = "%s/%d__%s/original_uploads/" % (directory, proj.id, proj.name )

        try:
            os.makedirs(files_path)
        except OSError:
            pass
        try:
            os.makedirs(uploads_path)
        except OSError:
            pass
        for cbh_file in CBHFlowFile.objects.filter(project=proj):
            #Make a new unique filename for the item
            if cbh_file.cbhcompoundmultiplebatch_set.count() == 0:
                try:
                    shutil.copy2(cbh_file.full_path, files_path + backup_filename(cbh_file.id, cbh_file.original_filename ))
                except IOError:
                    print "No such file %s" % cbh_file.full_path
            else:
                try:
                    shutil.copy2(cbh_file.full_path, uploads_path + backup_filename(cbh_file.id, cbh_file.original_filename ))
                except IOError:
                    print "No such file %s" % cbh_file.full_path

def backup_permissions(projects_list, directory):
    from cbh_core_api.resources import UserResource
    userres = UserResource()
    request_factory = RequestFactory()
    user = get_user_model().objects.filter(is_superuser=True)[0]
    for proj in projects_list:
        from django.contrib.auth.models import User, Permission
        from django.db.models import Q
        from cbh_core_model.models import PROJECT_PERMISSIONS, PERMISSION_CODENAME_SEPARATOR
        for perm_type, irrel, irel2 in PROJECT_PERMISSIONS:
            try:
                perm = Permission.objects.get(codename="%d%s%s" % (proj.id, PERMISSION_CODENAME_SEPARATOR, perm_type) )
            except Permission.DoesNotExist:
                continue
            users = User.objects.filter(Q(groups__permissions=perm) | Q(user_permissions=perm) ).distinct()
            project_path = "%s/%d__%s/" % (directory, proj.id, proj.name )
            try:
                os.makedirs(project_path)
            except OSError:
                pass
            newrequest = request_factory.get("/")
            newrequest.user = user
            newrequest.GET = newrequest.GET.copy()
            newrequest.GET["format"] = format
            newrequest.GET["offset"] = 0
            newrequest.GET["limit"] = 100000
            newrequest.GET["id__in"] = ",".join([str(u.id) for u in users])

            json_content = userres.get_list(newrequest).content
            path_to_write_json = project_path + "%s.json" % perm_type
            with open(path_to_write_json, "w") as text_file:
                text_file.write(json_content)





def delete_projects(projects_list):
    for proj in projects_list:
        from django.contrib.auth.models import User, Permission
        from cbh_core_model.models import PROJECT_PERMISSIONS, PERMISSION_CODENAME_SEPARATOR, CustomFieldConfig
        index_name = elasticsearch_client.get_project_index_name(proj.id)
        elasticsearch_client.delete_index(index_name)

        for perm_type, irrel, irel2 in PROJECT_PERMISSIONS:
            try:
                perm = Permission.objects.get(codename="%d%s%s" % (proj.id, PERMISSION_CODENAME_SEPARATOR, perm_type) )
                perm.delete()
            except :
                pass
        proj.custom_field_config.delete()
        print proj.name
        proj.delete()