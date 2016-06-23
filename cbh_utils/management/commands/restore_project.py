"""Give a project.json file, restore that project to the next available project ID"""
import json
from django.test import RequestFactory
import sys, os
from django.core.management.base import BaseCommand, CommandError
from django.contrib.auth import get_user_model
from django.db import transaction
import shortuuid

def restore_project_from_json(json_file):
    from cbh_core_api.resources import ChemregProjectResource, ProjectTypeResource
    with open(json_file) as f:
        data = f.read()
        jsondata = json.loads(data)
        old_id = jsondata["id"]
        del jsondata["id"]
        del jsondata["resource_uri"]
        del jsondata["created_by"]
        del jsondata["project_type"]["id"]
        del jsondata["project_type"]["resource_uri"]
        del jsondata["project_type"]["custom_field_config_template_id"]
        del jsondata["custom_field_config"]["created_by"]
        del jsondata["custom_field_config"]["id"]
        del jsondata["custom_field_config"]["resource_uri"]
        for field in jsondata["custom_field_config"]["project_data_fields"]:
            del field["id"]
            del field["resource_uri"]

        request_factory = RequestFactory()
        user = get_user_model().objects.filter(is_superuser=True)[0]
        request = request_factory.get("/")
        request.user = user
        request.GET = request.GET.copy()
        request.GET["name"] = jsondata["project_type"]["name"]
        resp = ProjectTypeResource().get_list(request)
        ptcont = json.loads(resp.content)
        if len(ptcont["objects"]) ==1:
            notsame = False
            project_type_to_check = ptcont["objects"][0]
            for fieldname in ["show_compounds", "saved_search_project_type", "plate_map_project_type"]:
                if project_type_to_check[fieldname] == jsondata["project_type"][fieldname]:
                    print ("checked %s and found same" % fieldname)
                else:
                    notsame = True
            if notsame:
                print ("project type found with same name that is not actually the same, quitting")
            #We found the project type is fine so just use the existing one
            jsondata["project_type"] = project_type_to_check["resource_uri"]
            #Otherwise just allow the script to create the given project type
        request = request_factory.post("/", content_type='application/json', data=json.dumps(jsondata))
        request.user = user
        resp = ChemregProjectResource().post_list(request)
        created_project = json.loads(resp.content)
        return (created_project, old_id)



def restore_permissions(new_project_id, directory):
    from cbh_core_model.models import PROJECT_PERMISSIONS
    from cbh_core_model.models import Project
    User = get_user_model()
    created_users = []
    for codename, stuff, stuf2 in PROJECT_PERMISSIONS:
        filename = "%s%s.json" % (directory, codename)
        with open(filename) as f:
            data = json.loads(f.read())
            for user_json in data["objects"]:
                u, created = User.objects.get_or_create(username=user_json["username"], defaults={"active": True})
                p = Project.objects.get(pk=new_project_id)
                p._add_instance_permissions_to_user_or_group(
            u,  codename)
                if created:
                    created_users.append(u)

    if len(created_users) > 0:
        print ("The following users have been created and will need to change password if using password login.")
        for cu in created_users:
            print (cu.username)



def get_compounds_json(directory):
    json_file = "%scompounds.json" % directory
    jsondata =  None
    with open(json_file) as f:
        data = f.read()
        jsondata = json.loads(data)
    return jsondata["objects"]



def restore_file_attachments(directory, compounds, project_json):
    from cbh_core_api.views import UploadView, CBHFlowFileResource
    attach_files = [filename for filename in os.listdir('%sattachments' % directory)]
    upload_url = project_json["flowjs_upload_url"]
    from django.conf import settings
    from importlib import import_module
    engine = import_module(settings.SESSION_ENGINE)
    
    
    for comp_json in compounds:
        for key, value in comp_json["custom_fields"].items():
            
            if isinstance(value, dict):
                if value.get("attachments", None):
                    for attachment in value["attachments"]:
                        url_bits = attachment["url"].split("/")
                        attach_id = url_bits[-1]
                        #pick out the filename
                        prefixed = [filename for filename in attach_files if filename.startswith(attach_id + "_")]
                        if len(prefixed) == 1:
                            with open('%sattachments/%s' % (directory,prefixed[0])) as f:
                                session_key = "None"
                                flowid = shortuuid.ShortUUID().random() + prefixed[0]
                                request_factory = RequestFactory()
                                user = get_user_model().objects.filter(is_superuser=True)[0]
                                request = request_factory.post(upload_url, {"file": f, "flowChunkNumber": 1, 
                                    "flowChunkSize": 22222222, 
                                    "flowCurrentChunkSize": 137227,
                                    "flowTotalSize": 137227,
                                    "flowFilename": prefixed[0],
                                    "flowIdentifier": flowid,
                                    "flowRelativePath": prefixed[0],
                                    "flowTotalChunks": 1})
                                request.session = engine.SessionStore(session_key)
                                uv = UploadView()
                                resp = uv.dispatch(request, project_id=project_json["id"])
                                request = request_factory.get("/")
                                request.session = engine.SessionStore(session_key)
                                cbhf = CBHFlowFileResource()
                                resp2 = cbhf.get_detail(request, identifier=flowid, pk=flowid)
                                jd = json.loads(resp2.content)
                                attachment["url"] = jd["download_uri"]
    return compounds

def restore_compounds(directory, compounds, project_json):
    from cbh_chem_api.resources import CBHCompoundBatchResource
    
    cbr = CBHCompoundBatchResource()
    request_factory = RequestFactory()
    user = get_user_model().objects.filter(is_superuser=True)[0]
    for comp_json in compounds:
        comp_json["warnings"]["original_uox_id"] = comp_json["uuid"]
        del comp_json["id"]
        del comp_json["resource_uri"]
        del comp_json["creator"]
        del comp_json["userfull"] 
        comp_json["project"] = {"pk": project_json["id"]}
        request = request_factory.post("/", data=json.dumps(comp_json), content_type='application/json',)
        print comp_json.keys()
        request.user = user
        resp = cbr.post_list(request)





class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('directory')



    def handle(self, *args, **options):
        with transaction.atomic():
            directory = options["directory"]
            project, old_id = restore_project_from_json("%sproject.json" % directory)
            restore_permissions(project["id"], directory)
            compounds = get_compounds_json(directory)
            compounds = restore_file_attachments(directory, compounds, project)
            restore_compounds(directory, compounds, project)




