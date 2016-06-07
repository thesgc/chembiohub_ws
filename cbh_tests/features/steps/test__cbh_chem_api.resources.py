from behave import given, when, then
import json
#from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType
from django.db import IntegrityError
#from django.contrib.auth.models import User


@given(u'I create a compound batch from a drawing as before')
def step_impl(context):
    context.execute_steps(u"""
       Given I create a project as before
        When I refresh the user object
        Given I have a compound batch with no structure
        and I add a valid molfile to my compound data and call it ctab
        and I add the project primary key to the compound data
        When I submit the compound by POST request
        Then a compound batch is created""")


@given(u'I save propane butane benzene and ethyl benzene via SMILES')
def step(context):
    context.execute_steps(u"""
        Given I set up the SMILES data
        When I validate propane butane benzene and ethyl benzene via SMILES
        Then the response from post validate list is accepted
        Then I request the multiple batch data until the response is OK
        then I can post to multi batch save and the response will eventually be created
        When I list compound batches in the system
        Then 4 compound batches have been created
    """)

@then(u'I request the multiple batch data until the response is OK')
def step_impl(context):
    status = 409
    while status != 200:
        from django.conf import settings
        url = "/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/get_part_processed_multiple_batch?current_batch=1" 
        resp = context.api_client.get(url, format='json')
        context.multibatch_response = resp
        status = resp.status_code
        context.test_case.assertTrue(status in [200, 409])
    context.test_case.assertTrue(status==200)

@then(u'I get a chemical search id in the response')
def step_impl(context):
    data = json.loads(context.chemical_search_response.content)
    context.test_case.assertTrue(isinstance(data["id"], basestring))
    context.test_case.assertTrue(len(data["id"]) > 10)
    context.chemical_search_id = data["id"]


@when(u'I run a chemical search')
def step(context):
    from django.conf import settings
    url = "/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches_v2/?chemical_search_id=" + context.chemical_search_id 
    print(url)
    resp = context.api_client.get(url, format='json')
    context.batches_response = resp

@then(u'I get 3 chemical search results')
def step(context):
    data = json.loads(context.batches_response.content)["objects"]
    print("found mols")
    print(len(data))
    context.test_case.assertTrue(len(data) == 3)

@then(u'I get 1 chemical search result')
def step(context):
    data = json.loads(context.batches_response.content)["objects"]
    print("found mols")
    print(len(data))
    context.test_case.assertTrue(len(data) == 1)




@when(u'I POST the molfile for propane as a flexmatch')
def step_impl(context):
    post_data = {"molfile":"Molecule from ChemDoodle Web Components\n\nhttp://www.ichemlabs.com\n  3  2  0  0  0  0            999 V2000\n   -0.8660   -0.2500    0.0000 C   0  0  0  0  0  0\n    0.0000    0.2500    0.0000 C   0  0  0  0  0  0\n    0.8660   -0.2500    0.0000 C   0  0  0  0  0  0\n  1  21  0     0  0\n  2  31  0     0  0\nM  END","query_type":"flexmatch"}
    from django.conf import settings
    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_chemical_search/", format='json', data=post_data)
    context.chemical_search_response = resp


@when(u'I POST the molfile for propane as a with_substructure')
def step_impl(context):
    post_data = {"molfile":"Molecule from ChemDoodle Web Components\n\nhttp://www.ichemlabs.com\n  3  2  0  0  0  0            999 V2000\n   -0.8660   -0.2500    0.0000 C   0  0  0  0  0  0\n    0.0000    0.2500    0.0000 C   0  0  0  0  0  0\n    0.8660   -0.2500    0.0000 C   0  0  0  0  0  0\n  1  21  0     0  0\n  2  31  0     0  0\nM  END","query_type":"with_substructure"}
    from django.conf import settings
    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_chemical_search/", format='json', data=post_data)
    context.chemical_search_response = resp






@when(u'I validate propane butane benzene and ethyl benzene via SMILES')
def step(context):
    from django.conf import settings
    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/validate_list/", format='json', data=context.post_data)
    context.val_response = resp

@given(u'I set up the SMILES data')
def step(context):
    context.execute_steps(u"""
        Given I create a project as before
        When I refresh the user object
        Given I have a compound batch with no structure
        and I add the project key to the compound data
        and I set the state to validate
    """)
    context.post_data["smilesdata"] = """CCc1ccccc1
CCC
CCCC
c1ccccc1"""

@when(u'I save the post validate multiple batch')
def step(context):
    from django.conf import settings
    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/multi_batch_save/", format='json', data=context.post_data)




@when(u'I request the compound batch with ID 1 from the get_detail api')
def step_impl(context):
    from django.conf import settings
    resp = context.api_client.get("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches_v2/1" , format='json')
    context.batch_response = resp

@then(u'the batch response is OK')
def step_impl(context):
    context.test_case.assertHttpOK(context.batch_response)


@when(u'I list compound batches in the system')
def step_impl(context):
    from django.conf import settings
    resp = context.api_client.get("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches_v2/", format='json')
    context.batches_response = resp

@then(u'the created compound batch has a uox id in the uuid field')
def step_impl(context):
    from django.conf import settings
    data = json.loads(context.batches_response.content)["objects"]
    print("final")
    print(data)
    context.test_case.assertEquals(len(data), 1)

    
    context.test_case.assertTrue(data[0]["uuid"].startswith(settings.ID_PREFIX))


@then("the created compound batch has a multiple_batch_id")
def step(context):
    data = json.loads(context.batches_response.content)["objects"]
    context.test_case.assertEquals(len(data), 1)
    context.test_case.assertEquals(data[0]["multiple_batch_id"], 1)




@then(u'the compound batch list response is OK')
def step_impl(context):
    context.test_case.assertHttpOK(context.batches_response)

@then(u'the saved search response is OK')
def step_impl(context):
    context.test_case.assertHttpOK(context.batches_response)


@then(u'I see no compound batches')
def step_impl(context):
    data = json.loads(context.batches_response.content)["objects"]
    context.test_case.assertEquals(len(data), 0)
    print (context.batches_response.content)



@when(u'I list saved searches in the system')
def step_impl(context):
    from django.conf import settings
    resp = context.api_client.get("/" + settings.WEBSERVICES_NAME + "/cbh_saved_search/", format='json')
    context.batches_response = resp


@then(u'I see my saved search')
def step_impl(context):
    data = json.loads(context.batches_response.content)["objects"]
    context.test_case.assertEquals(len(data), 1)
    context.test_case.assertEquals(data[0]["custom_fields"]["alias"], "My new saved search")





@then(u'The saved search index response is OK')
def step_impl(context):
    context.test_case.assertHttpOK(context.saved_search_index_resp)


@given(u'A URL to redirect the user to and a GET request URL and an alias and description for my saved search')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.saved_search_data = {
        "project": {"pk" : project_json["id"]},
        "custom_fields" : {"alias": "My new saved search",
                            "url" : "http://localhost:9000/search"},
        "uncurated_fields":{},
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


@given(u'I add the project key')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.saved_search_data["project_key"] =  project_json["project_key"]



@given("I create a saved search as before")
def step_impl(context):
    context.execute_steps(u"""
        Given I have loaded the fixtures for project types
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 4
        Given I create a project JSON by adding the saved search project type and the static custom field config schema for saved searches
        When I POST a project to cbh_projects
        Then the project is created
        When I refresh the user object
        Given A URL to redirect the user to and a GET request URL and an alias and description for my saved search
        When I send the search by POST request
        Then The saved search response is created

        """)

        
@given(u'I have a compound batch with no structure')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.post_data = {
        "project": project_json["resource_uri"],
        "custom_fields" : {"test": "My test field"
                         },
        "uncurated_fields":{},
        "warnings" :{}, 
        "properties" :{}, 
        "errors" :{}
    }

@given(u'I set the restricted field to the value foo')
def step(context):
    context.post_data["restricted"] = "foo"





@then(u'the response from post validate list is accepted')
def step_impl(context):
    context.test_case.assertHttpAccepted(context.val_response )
    print (context.val_response)
    context.valdata = json.loads(context.val_response.content)


@then(u'4 compound batches have been created')
def step(context):
    data = json.loads(context.batches_response.content)["objects"]
    context.test_case.assertEquals(len(data), 4)




@given(u"A compound batch is created from a drawing as before")
def step_impl(context):
    context.execute_steps(u"""
        Given I create a project as before
        When I refresh the user object
        Given I have a compound batch with no structure
        and I add a valid molfile to my compound data and call it ctab
        and I add the project key to the compound data
        and I set the state to validate
        When I submit the compound to POST validate drawn
        then the response from post validate drawn is accepted
        Then I request the multiple batch data until the response is OK 

        then I can post to multi batch save and the response will eventually be created
        """)


@then(u'I can post to multi batch save and the response will eventually be created')
def step_impl(context):
    status = 202

    while status == 202:
        from django.conf import settings
        data_to_post = json.loads(context.multibatch_response.content)
        context.multibatch_response = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/multi_batch_save/", format='json', data= data_to_post)
        import time
        status = context.multibatch_response.status_code
        time.sleep(1)

    context.test_case.assertHttpCreated(context.multibatch_response)



@when(u'a compound exists in the project with this restricted field')
def step(context):
    context.execute_steps(u"""
        Given I have a compound batch with no structure
        Given I set the restricted field to the value foo
        Given I add the project primary key to the compound data
        When I submit the compound by POST request
        Then a compound batch is created        
        """)


@when(u'I submit the compound by POST request')
def step_impl(context):
    from django.conf import settings
    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches_v2/", format='json', data=context.post_data)
    context.compound_response = resp
    print (resp.content)


@then(u'the output data from search does not contain the restricted field')
def step(context):
    data = json.loads(context.batches_response.content)
    batch = data["objects"][0]
    context.test_case.assertTrue("restricted" not in batch["custom_fields"])
    context.test_case.assertTrue("test" in batch["custom_fields"])


@then(u'a compound batch is created')
def step_impl(context):
    context.test_case.assertHttpCreated(context.compound_response)


#The below scenarios are inconsistent and need an update
@given(u'I add the project key to the compound data')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.post_data["project_key"] =  project_json["project_key"]


@given(u'I add the project primary key to the compound data')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.post_data["project"] =  {"pk" : project_json["id"] }



@given(u'I set the state to validate')
def step_impl(context):
    project_json = json.loads(context.project_response.content)
    context.post_data["validate"] =  "validate"





@then(u'the blinded batch ID is generated')
def step_impl(context):
    from django.conf import settings
    comp = json.loads(context.compound_response.content)
    print (comp["blinded_batch_id"])
    context.test_case.assertTrue(comp["blinded_batch_id"].startswith(settings.ID_PREFIX))



@given('I add a valid molfile to my compound data and call it ctab')
def step(context):
    context.post_data["ctab"] = """


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
    "project_key": project_key,
    "struccol":"",
    "state":"validate"}
    from django.conf import settings
    context.val_response = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/cbh_compound_batches/validate_files/", 
        format='json', 
        data=validate_file_data)
