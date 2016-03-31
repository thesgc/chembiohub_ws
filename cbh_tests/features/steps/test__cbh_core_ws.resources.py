from behave import given, when, then
import json
from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType
from django.db import IntegrityError
from django.contrib.auth.models import User


@then("I can list the users on the system")
def users(context):
    from django.conf import settings
    resp = context.api_client.client.get("/" +settings.WEBSERVICES_NAME + "/users/")
    context.test_case.assertHttpOK(resp)

@then("I can list the projects types on the system and there are 4")
def project_types(context):
    from django.conf import settings
    resp = context.api_client.client.get("/" + settings.WEBSERVICES_NAME + "/cbh_project_types/")
    context.test_case.assertHttpOK(resp)
    json_content = json.loads(resp.content)
    context.test_case.assertEqual(len(json_content["objects"]), 4)
    context.project_types = json_content["objects"]


@then("I can list the data types on the system and there are 4")
def data_types(context):
    from django.conf import settings
    resp = context.api_client.client.get("/" + settings.WEBSERVICES_NAME + "/datastore/cbh_data_types/")
    context.test_case.assertHttpOK(resp)
    json_content = json.loads(resp.content)
    context.test_case.assertEqual(len(json_content["objects"]), 4)
    context.data_types = json_content["objects"]


@given("testuser")
def step(context):
    pass


@when("I log in testuser")
def logintestuser(context):
    context.api_client.client.login(username="testuser", password="testuser")

@given("I log out")
def logout(context):
    context.api_client.client.logout()

@then("I can see the skinning configuration")
def skinning(context):
    from django.conf import settings
    resp = context.api_client.get("/" +settings.WEBSERVICES_NAME + "/cbh_skinning/")
    context.test_case.assertHttpOK(resp)

    





@then(u'the default project type has an id of 3')
def step_impl(context):
    id = -1
    for ptype in context.project_types:
        if ptype["set_as_default"]:
            id=ptype["id"]
    context.test_case.assertEqual(id, 3)


@when(u'I take the project type model with id 1 and set it to default')
def step_impl(context):
    ptype = ProjectType.objects.get(id=1)
    ptype.set_as_default = True
    ptype.save()


@then(u'the previous default project type has been removed by the system and only the new one with id 1 remains')
def step_impl(context):
    qs = ProjectType.objects.filter(set_as_default=True)
    context.test_case.assertEqual(qs.count(), 1)
    context.test_case.assertEqual(qs[0].id, 1)


@then(u'the project type with id 4 has the saved_search_project_type')
def step_impl(context):
    id = -1
    for ptype in context.project_types:
        if ptype["saved_search_project_type"]:
            id=ptype["id"]
    context.test_case.assertEqual(id, 4)

@then(u'the project type with id 4 has a static project template to allow saved searches to be saved as compound batches')
def step_impl(context):
    project_template = None
    for ptype in context.project_types:
        if ptype["saved_search_project_type"]:
            project_template = json.dumps(ptype["project_template"]["custom_field_config"]["project_data_fields"])
    #This variable is a copy of the variable in ProjectType.SAVED_SEARCH_TEMPLATE, if you reimplement you must also change the test as an insurance policy for the front end code
    STATIC_SHOUlD_BE_UNLESS_SYSTEM_CHANGED = [{
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
    context.test_case.assertEqual(project_template, json.dumps(STATIC_SHOUlD_BE_UNLESS_SYSTEM_CHANGED))












