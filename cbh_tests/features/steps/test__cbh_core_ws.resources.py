from behave import given, when, then
import json
from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType
from cbh_datastore_model.models import DataPoint, DataPointClassification, DataPointClassificationPermission
from django.db import IntegrityError
from django.contrib.auth.models import User


@then("I can list the users on the system")
def users(context):
    from django.conf import settings
    resp = context.api_client.client.get("/" +settings.WEBSERVICES_NAME + "/users/")
    context.test_case.assertHttpOK(resp)

@then("I can list the projects types on the system and there are 3")
def project_types(context):
    from django.conf import settings
    resp = context.api_client.client.get("/" + settings.WEBSERVICES_NAME + "/cbh_project_types/")
    context.test_case.assertHttpOK(resp)
    json_content = json.loads(resp.content)
    context.test_case.assertEqual(len(json_content["objects"]), 3)
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

    
