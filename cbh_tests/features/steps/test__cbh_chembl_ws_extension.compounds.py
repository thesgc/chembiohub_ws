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
        "project": 
        "customFields" : {"alias": "My new saved search",
                            "url" : "http://localhost:9000/search"}
    }

@when(u'I send the search by POST request')
def step_impl(context):
    from django.conf import settings

    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/post_saved_search/", format='json', data=context.saved_search_data)
    context.saved_search_response = resp
    print(resp.content)

    

@then(u'The saved search response is created')
def step_impl(context):
    raise NotImplementedError(u'STEP: Then The saved search response is created')

@then(u'the creator of the saved search is set as its owner')
def step_impl(context):
    raise NotImplementedError(u'STEP: Then the creator of the saved search is set as its owner')
