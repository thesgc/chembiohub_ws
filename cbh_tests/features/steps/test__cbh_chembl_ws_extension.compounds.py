from behave import given, when, then
import json
from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType
from cbh_datastore_model.models import DataPoint, DataPointClassification, DataPointClassificationPermission
from django.db import IntegrityError
from django.contrib.auth.models import User








@given(u'A URL to redirect the user to and a GET request URL and an alias and description for my saved search')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.saved_search_data = {
        "project": project_json["resource_uri"],
        "customFields" : {"alias": "My new saved search",
                            "url" : "http://localhost:9000/search"},

    }

@when(u'I send the search by POST request')
def step_impl(context):
    from django.conf import settings

    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_saved_search/", format='json', data=context.saved_search_data)
    context.saved_search_response = resp
    print(resp.content)



@then(u'The saved search response is created')
def proj_create(context):
    context.test_case.assertHttpCreated(context.saved_search_response)


@given(u'I add the blinded batch id as EMPTY_STRING')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.saved_search_data["blinded_batch_id"] =  "EMPTY_STRING"


@given(u'I add the project key')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.saved_search_data["project_key"] =  project_json["project_key"]



@given("I create a saved search as before")
def step_impl(context):
    context.execute_steps(u"""
        Given I have loaded the fixtures for project types and data types
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 4
        Given I create a project JSON by adding the saved search project type and the static custom field config schema for saved searches
        When I POST a project to cbh_projects
        Then the project is created
        Given A URL to redirect the user to and a GET request URL and an alias and description for my saved search
        and I add the blinded batch id as EMPTY_STRING
        and I add the project key
        When I send the search by POST request
        Then The saved search response is created
        """)

        