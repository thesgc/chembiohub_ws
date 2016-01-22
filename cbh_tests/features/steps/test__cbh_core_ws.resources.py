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

    
