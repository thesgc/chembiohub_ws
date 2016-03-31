from behave import given, when, then
import json
from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType
from django.db import IntegrityError
from django.contrib.auth.models import User













@given(u'I create a compound batch from a drawing as before')
def step_impl(context):
    context.execute_steps(u"""
        Given I create a project as before
        When I refresh the user object
        Given I have a compound batch with no structure
        and I add a valid molfile to my compound data and call it sketch
        and I add the project key to the compound data
        and I set the type of the request to sketch
        and I set the state to validate
        When I submit the compound to POST validate drawn
        then the response from post validate drawn is accepted
        when I take the response from post validate drawn and post it to multi batch save
        then the response from multi batch save is created""")


@when(u'I request the compound batch with ID 1 from the get_detail api')
def step_impl(context):
    from django.conf import settings
    resp = context.api_client.get("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/1" , format='json')
    context.batch_response = resp

@then(u'the batch response is OK')
def step_impl(context):
    context.test_case.assertHttpOK(context.batch_response)


@when(u'I list compound batches in the system with get_list_elasticsearch')
def step_impl(context):
    from django.conf import settings
    resp = context.api_client.get("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/get_list_elasticsearch/", format='json')
    context.batches_response = resp

@then(u'the created compound batch has a uox id in the chemblId field')
def step_impl(context):
    from django.conf import settings
    data = json.loads(context.batches_response.content)["objects"]
    context.test_case.assertEquals(len(data), 1)
    print(data)
    context.test_case.assertTrue(data[0]["chemblId"].startswith(settings.ID_PREFIX))


@then("the created compound batch has a multipleBatchId")
def step(context):
    data = json.loads(context.batches_response.content)["objects"]
    context.test_case.assertEquals(len(data), 1)
    context.test_case.assertEquals(data[0]["multipleBatchId"], 1)




@then(u'the get list elasticsearch response is OK')
def step_impl(context):
    context.test_case.assertHttpOK(context.batches_response)

@then(u'I see no compound batches')
def step_impl(context):
    data = json.loads(context.batches_response.content)["objects"]
    context.test_case.assertEquals(len(data), 0)
    print (context.batches_response.content)



@when(u'I list saved searches in the system with get_list_elasticsearch')
def step_impl(context):
    from django.conf import settings
    resp = context.api_client.get("/" + settings.WEBSERVICES_NAME + "/cbh_saved_search/get_list_elasticsearch/", format='json')
    context.batches_response = resp


@then(u'I see my saved search')
def step_impl(context):
    data = json.loads(context.batches_response.content)["objects"]
    context.test_case.assertEquals(len(data), 1)
    context.test_case.assertEquals(data[0]["customFields"]["alias"], "My new saved search")



@when(u'I reindex the saved search')
def step_impl(context):
    from django.conf import settings
    saved_search_json = json.loads(context.saved_search_response.content)
    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_saved_search/reindex_compound/", format='json', data={"id": saved_search_json["id"]})
    context.saved_search_index_resp = resp


@then(u'The saved search index response is OK')
def step_impl(context):
    context.test_case.assertHttpOK(context.saved_search_index_resp)


@given(u'A URL to redirect the user to and a GET request URL and an alias and description for my saved search')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.saved_search_data = {
        "project": project_json["resource_uri"],
        "customFields" : {"alias": "My new saved search",
                            "url" : "http://localhost:9000/search"},
        "uncuratedFields":{},
        "warnings" :{}, "properties" :{},  "errors" :{}
    }


@when(u'I send the search by POST request')
def step_impl(context):
    from django.conf import settings
    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_saved_search/", format='json', data=context.saved_search_data)
    context.saved_search_response = resp



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
        When I reindex the saved search 
        Then The saved search index response is OK

        """)

        
@given(u'I have a compound batch with no structure')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.post_data = {
        "project": project_json["resource_uri"],
        "customFields" : {"test": "My test field"
                         },
        "uncuratedFields":{},
        "warnings" :{}, 
        "properties" :{}, 
        "errors" :{}
    }


@given(u'I set the type of the request to sketch')
def step(context):
    context.post_data["type"] = "sketch"



@given(u'I add the blinded batch id to my compound POST data as EMPTY_STRING')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.post_data["blinded_batch_id"] =  "EMPTY_STRING"




@when(u'I submit the compound to POST validate drawn')
def step_impl(context):
    from django.conf import settings
    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/validate_drawn/", format='json', data=context.post_data)
    context.val_response = resp

@then(u'the response from post validate drawn is accepted')
def step_impl(context):
    context.test_case.assertHttpAccepted(context.val_response )
    print (context.val_response)
    context.valdata = json.loads(context.val_response.content)





@when(u'I take the response from post validate drawn and post it to multi batch save')
def step_impl(context):
    from django.conf import settings
    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/multi_batch_save/", format='json', data=context.valdata )
    context.multibatch_response = resp


@given(u"A compound batch is created from a drawing as before")
def step_impl(context):
    context.execute_steps(u"""
        Given I create a project as before
        When I refresh the user object
        Given I have a compound batch with no structure
        and I add a valid molfile to my compound data and call it sketch
        and I add the project key to the compound data
        and I set the type of the request to sketch
        and I set the state to validate
        When I submit the compound to POST validate drawn
        then the response from post validate drawn is accepted
        when I take the response from post validate drawn and post it to multi batch save
        then the response from multi batch save is created
        """)


@then(u'the response from multi batch save is created')
def step_impl(context):
    context.test_case.assertHttpCreated(context.multibatch_response)
    print (context.multibatch_response.content)



@when(u'I submit the compound by POST request')
def step_impl(context):
    from django.conf import settings
    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/", format='json', data=context.post_data)
    context.compound_response = resp


@then(u'a compound batch is created')
def step_impl(context):
    context.test_case.assertHttpCreated(context.compound_response)


@given(u'I add the project key to the compound data')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.post_data["project_key"] =  project_json["project_key"]


@given(u'I set the state to validate')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.post_data["validate"] =  "validate"





@when(u'I reindex the compound')
def step_impl(context):
    from django.conf import settings
    comp = json.loads(context.compound_response.content)
    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_saved_search/reindex_compound/", format='json', data={"id": comp["id"]})
    context.compound_index_response = resp


@then(u'the blinded batch ID is generated')
def step_impl(context):
    from django.conf import settings
    comp = json.loads(context.compound_response.content)
    context.test_case.assertTrue(comp["blindedBatchId"].startswith(settings.ID_PREFIX))



@given('I add a valid molfile to my compound data and call it sketch')
def step(context):
    context.post_data["sketch"] = """


  8  8  0  0  0  0            999 V2000
    0.0000    1.0000    0.0000 C   0  0  0  0  0  0
    0.8660    0.5000    0.0000 C   0  0  0  0  0  0
    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0
    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0
   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0
   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0
    0.0000    2.0000    0.0000 O   0  0  0  0  0  0
    0.0000   -2.0000    0.0000 O   0  0  0  0  0  0
  1  21  0     0  0
  2  32  0     0  0
  3  41  0     0  0
  4  51  0     0  0
  5  62  0     0  0
  6  11  0     0  0
  1  72  0     0  0
  4  82  0     0  0
M  END"""






@then(u'the response from post validate files is accepted')
def step_impl(context):
    context.test_case.assertHttpAccepted(context.val_response )
    print (context.val_response)
    context.valdata = json.loads(context.val_response.content)








#File upload stuff


@when("I validate the compounds file")
def step(context):
    project_key = json.loads(context.project_response.content )["project_key"]
    validate_file_data = {"file_name": context.current_file_name,
    "multiplebatch":None,
    "type":"file",
    "fileextension":context.current_file_extension,
    "projectKey": project_key,
    "struccol":"",
    "state":"validate"}
    from django.conf import settings
    context.val_response = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/validate_files/", 
        format='json', 
        data=validate_file_data)
