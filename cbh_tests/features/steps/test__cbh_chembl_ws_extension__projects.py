from behave import given, when, then
import json
from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType
from cbh_datastore_model.models import DataPoint, DataPointClassification, DataPointClassificationPermission
from django.db import IntegrityError
from django.contrib.auth.models import User


@then("I can list the projects on the system")
def projects(context):
    from django.conf import settings
    resp = context.api_client.client.get("/" + settings.WEBSERVICES_NAME + "/cbh_projects/")
    context.test_case.assertHttpOK(resp)







@when("I POST a project to cbh_projects")
def project_create(context):
    from django.conf import settings


    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_projects/", format='json', data=context.project_json)
    context.project_response = resp


@then("I can list the projects types on the system and there are 3")
def project_types(context):
    from django.conf import settings
    resp = context.api_client.client.get("/" + settings.WEBSERVICES_NAME + "/cbh_project_types/")
    context.test_case.assertHttpOK(resp)
    json_content = json.loads(resp.content)
    context.test_case.assertEqual(len(json_content["objects"]), 3)
    context.project_types = json_content["objects"]


@given("I create a project JSON by adding one of these project types and some custom fields and a name")
def build_project_json(context):
    context.project_json = {
        "project_type": context.project_types[0],
            "custom_field_config": {
                "project_data_fields": [{
                    "required": False,
                    "field_type": "char",
                    "open_or_restricted": "open",
                    "name": "Test Field"
                }],
                "name": "Test Project 2__config"
            },
            "name": "Test Project 2"
        }



@then("The project is created")
def proj_create(context):
    context.test_case.assertHttpCreated(context.project_response)


@then("The project is not created and there is a conflict")
def proj_create(context):
    context.test_case.assertHttpConflict(context.project_response)




@given("testuser has the cbh_core_model.add_project permission")
def add_proj_perms(context):
    from django.contrib.auth.models import Permission

    perm = Permission.objects.get_by_natural_key(codename="add_project", app_label="cbh_core_model", model="project")
    context.user.user_permissions.add(perm)
    context.user.save()