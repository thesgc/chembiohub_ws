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
    context.projects_on_system = json.loads(resp.content)["objects"]


@when("I patch the updated first project in the list back to the system")
def project_patch(context):
    from django.conf import settings
    resp = context.api_client.patch(context.projects_on_system[0]["resource_uri"], data=context.projects_on_system[0])
    context.updated_project_response = resp



@then("project update response is accepted")
def accepted(context):
    context.test_case.assertHttpAccepted(context.updated_project_response)



@then("the name has changed to Foo and there is one project in the list")
def name_changed(context):
    context.test_case.assertEqual(len(context.projects_on_system), 1)
    context.test_case.assertEqual(context.projects_on_system[0]["name"], "Foo")


@given("I take the first project in the list and change the name to be equal to the second one")
def project_name_change_for_conflict(context):
    context.projects_on_system[0]["name"] = context.projects_on_system[1]["name"]



@given("I take the first project in the list and change the name to Foo")
def project_name_change(context):
    context.projects_on_system[0]["name"] = "Foo"


@when("I POST a project to cbh_projects")
def project_create(context):
    from django.conf import settings


    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_projects/", format='json', data=context.project_json)
    context.project_response = resp
    print(resp.content)


@then("I can list the projects types on the system and there are 3")
def project_types(context):
    from django.conf import settings
    resp = context.api_client.client.get("/" + settings.WEBSERVICES_NAME + "/cbh_project_types/")
    context.test_case.assertHttpOK(resp)
    json_content = json.loads(resp.content)
    context.test_case.assertEqual(len(json_content["objects"]), 3)
    context.project_types = json_content["objects"]


@given("I create a project JSON by adding one of these project types and some custom fields and Bar as a name")
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
                "name": "Bar__config"
            },
            "name": "Bar"
        }
    print (context.project_json )



@given("I create a project JSON by adding one of these project types and some custom fields and Foo as a name")
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
                "name": "Foo_config"
            },
            "name": "Foo"
        }
    print (context.project_json )




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

@then("the project update response is unauthorized")
def proj_update_unauthorized(context):
    context.test_case.assertHttpUnauthorized(context.updated_project_response)

@then("the project is created not as unauthorized")
def proj_unauthorized(context):
    context.test_case.assertHttpUnauthorized(context.project_response)


@then("the project update response is conflict")
def proj_conflict(context):
    context.test_case.assertHttpConflict(context.updated_project_response)
