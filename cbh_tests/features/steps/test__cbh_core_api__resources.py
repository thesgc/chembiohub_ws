from behave import given, when, then
import json
from django.db import IntegrityError
import json


@then("all required fields are in the skinning configuration")
def step(context):
    for field in ["file_errors_from_backend","data_manager_email","chem_query_schemaform","id","enable_smiles_input","tabular_data_schema","hide_schemaform","field_type_choices","query_schemaform","filters_objects","savedsearch_schemaform","data_manager_name","sorts_applied","result_alias","project_alias","hides_applied","instance_alias","filters_applied","hide_objects","sort_objects","sort_schemaform","resource_uri"]:
        context.test_case.assertTrue(field in context.skinning["objects"][0])



@then("I can list the users on the system")
def users(context):
    from django.conf import settings
    context.userresp = context.api_client.client.get("/" +settings.WEBSERVICES_NAME + "/users/")
    context.test_case.assertHttpOK(context.userresp)

@then("I can list the projects types on the system and there are 4")
def project_types(context):
    from django.conf import settings
    resp = context.api_client.client.get("/" + settings.WEBSERVICES_NAME + "/cbh_project_types/")
    context.test_case.assertHttpOK(resp)
    json_content = json.loads(resp.content)
    context.test_case.assertEqual(len(json_content["objects"]), 4)
    context.project_types = json_content["objects"]


@then("the tabular schema included in fields all exist in the schema")
def step(context):
    counter = 0
    for key, value in context.tabular_schema["included_in_tables"].items():
        for fieldname in value["default"]:
            print (fieldname)
            context.test_case.assertTrue(fieldname in context.tabular_schema["schema"])
            counter += 1
    context.test_case.assertTrue(counter > 0)

@then("all the tabular fields have a knownBy value")
def step(context):
    counter = 0
    for key, value in context.tabular_schema["included_in_tables"].items():
        for fieldname in value["default"]:
            print (fieldname)
            knownBy = context.tabular_schema["schema"][fieldname]["knownBy"]
            context.test_case.assertTrue(len(knownBy) > 0)
            context.test_case.assertTrue(isinstance(knownBy, basestring))
            counter += 1
    context.test_case.assertTrue(counter > 0)

@then("all the tabular Data Schema included in fields have a project specific schema with a string ID inside it if they are custom fields")
def step(context):
    counter = 0
    for key, value in context.tabular_schema["included_in_tables"].items():
        for fieldname in value["default"]:
            if fieldname.startswith("custom_fields") and fieldname != "custom_fields.wells":
                proj_spec = context.tabular_schema["schema"][fieldname]["project_specific_schema"]
                context.test_case.assertTrue("1" in proj_spec)
                counter += 1
    context.test_case.assertTrue(counter > 0)




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
    context.skinning = json.loads(resp.content)
    context.test_case.assertHttpOK(resp)




@then(u'the default project type has an id of 2')
def step_impl(context):
    id = -1
    for ptype in context.project_types:
        if ptype["set_as_default"]:
            id=ptype["id"]
    context.test_case.assertEqual(id, 2)


@when(u'I take the project type model with id 1 and set it to default')
def step_impl(context):
    from cbh_core_model.models import ProjectType
    ptype = ProjectType.objects.get(id=1)
    ptype.set_as_default = True
    ptype.save()


@then(u'the previous default project type has been removed by the system and only the new one with id 1 remains')
def step_impl(context):
    from cbh_core_model.models import ProjectType
    qs = ProjectType.objects.filter(set_as_default=True)
    context.test_case.assertEqual(qs.count(), 1)
    context.test_case.assertEqual(qs[0].id, 1)


@then(u'the project type with id 3 has the saved_search_project_type')
def step_impl(context):
    id = -1
    for ptype in context.project_types:
        if ptype["saved_search_project_type"]:
            id=ptype["id"]
    context.test_case.assertEqual(id, 3)

@then(u'the project type with id 3 has a static project template to allow saved searches to be saved as compound batches')
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









@then("the first project in the list contains the restricted field in the schemata")
def step(context):
    fields = context.projects_on_system[0]["custom_field_config"]["project_data_fields"]
    has_restricted = False
    for field in fields:
        if field["name"] == "restricted":
            has_restricted = True
    context.test_case.assertTrue(has_restricted)

    has_restricted = False
    for field_path in context.tabular_schema["schema"]:
        if field_path == "custom_fields.restricted":
            has_restricted = True

    context.test_case.assertTrue(has_restricted)


@then("the first project in the list does not contain the restricted field in the schemata")
def step(context):
    fields = context.projects_on_system[0]["custom_field_config"]["project_data_fields"]
    has_restricted = False
    for field in fields:
        if field["name"] == "restricted":
            has_restricted = True
    context.test_case.assertFalse(has_restricted)

    has_restricted = False
    for field_path in context.tabular_schema["schema"]:
        if field_path == "custom_fields.restricted":
            has_restricted = True

    context.test_case.assertFalse(has_restricted)



@when("I request the owner permissions for the first project in the list")
def step(context):
    from django.conf import settings
    from cbh_core_model.models import  PERMISSION_CODENAME_SEPARATOR
    codename = "%d%sowner" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )
    context.permissions_response = context.api_client.client.get("/" + settings.WEBSERVICES_NAME + "/cbh_permissions/" + codename)  
    

@then("the response permission response is OK")
def step(context):
    context.test_case.assertHttpOK(context.permissions_response)

@then("the permissions response is unauthorized")
def step(context):
    context.test_case.assertHttpUnauthorized(context.permissions_response)


@when("I request the editor permissions for the first project in the list")
def step(context):
    from django.conf import settings
    from cbh_core_model.models import  PERMISSION_CODENAME_SEPARATOR
    codename = "%d%seditor" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )
    context.permissions_response = context.api_client.client.get("/" + settings.WEBSERVICES_NAME + "/cbh_permissions/" + codename)  
    


@when(u'I patch back the editor permissions with anotheruser added')
def step_impl(context):
    context.execute_steps(u"""
        Then I can list the users on the system
        """)
    user_list = json.loads(context.userresp.content)["objects"]
    
    for user in user_list:
        if user["username"] == "anotheruser":
            uri = user["resource_uri"]

    perms = json.loads(context.permissions_response.content)
    perms["users"].append(uri)

    from django.conf import settings
    from cbh_core_model.models import  PERMISSION_CODENAME_SEPARATOR
    codename = "%d%seditor" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )
    context.permissions_response = context.api_client.patch("/" + settings.WEBSERVICES_NAME + "/cbh_permissions/" + codename, data=perms)  





@when(u'I patch back the owner permissions with anotheruser added')
def step_impl(context):
    context.execute_steps(u"""
        Then I can list the users on the system
        """)
    user_list = json.loads(context.userresp.content)["objects"]
    
    for user in user_list:
        if user["username"] == "anotheruser":
            uri = user["resource_uri"]

    perms = json.loads(context.permissions_response.content)
    perms["users"].append(uri)

    from django.conf import settings
    from cbh_core_model.models import  PERMISSION_CODENAME_SEPARATOR
    codename = "%d%sowner" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )
    context.permissions_response = context.api_client.patch("/" + settings.WEBSERVICES_NAME + "/cbh_permissions/" + codename, data=perms)  


@then(u'the permissions response is accepted')
def step_impl(context):
    context.test_case.assertHttpAccepted(context.permissions_response)


@then("the permissions response shows that I am an owner")
def step(context):
    context.execute_steps(u"""
        Then I can list the users on the system
        """)
    user_list = json.loads(context.userresp.content)["objects"]

    for user in user_list:
        if user["username"] == "testuser":
            uri = user["resource_uri"]

    perms = json.loads(context.permissions_response.content)

    #The uri for testuser should be equal to the list of owners
    context.test_case.assertEqual(perms["users"], [uri,])

@when("I patch back the owner permissions as empty")
def step(context):
    perms = json.loads(context.permissions_response.content)
    perms["users"] = []
    from django.conf import settings
    from cbh_core_model.models import  PERMISSION_CODENAME_SEPARATOR
    print(perms)
    codename = "%d%sowner" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )
    context.permissions_response = context.api_client.patch("/" + settings.WEBSERVICES_NAME + "/cbh_permissions/" + codename, 
                                                                   format='json',
                                                                   data=perms)




@then("I can list the projects on the system")
def projects(context):
    from django.conf import settings
    resp = context.api_client.client.get("/" + settings.WEBSERVICES_NAME + "/cbh_projects/")
    context.test_case.assertHttpOK(resp)
    context.projects_on_system = json.loads(resp.content)["objects"]#
    context.tabular_schema = json.loads(resp.content)["tabular_data_schema"]


@then("the upload URL from the first project in the list points to the right place")
def test_url(context):
    from django.conf import settings
    context.upload_url = context.projects_on_system[0]["flowjs_upload_url"]
    context.test_case.assertEqual(context.upload_url, "/" + settings.WEBSERVICES_NAME + "/flowv2/upload/1/")



@when("I patch the updated first project in the list back to the system")
def project_patch(context):
    resp = context.api_client.patch(context.projects_on_system[0]["resource_uri"], format='json', data=context.projects_on_system[0])
    context.updated_project_response = resp

@given("I add a restricted field")
def step(context):
    new_field = { 
                    "name" : "restricted", 
                    "open_or_restricted" : "restricted",  
                    "field_type" : "char", 
                    "description" : "" 
                }
    context.projects_on_system[0]["custom_field_config"]["project_data_fields"].append(new_field)


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


@given("I create a project as before")
def step(context):
    context.execute_steps(u"""
        Given I have loaded the fixtures for project types
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 4
        Given I create a project JSON by adding one of these project types and some custom fields and Bar as a name
        When I POST a project to cbh_projects
        Then the project is created
        Then I can list the projects on the system
        """)




@given("I create a project JSON by adding the saved search project type and the static custom field config schema for saved searches")
def build_project_json(context):
    for pType in context.project_types:
        if pType["name"] == "Saved Search":
            context.project_json = {
                "project_type": pType,
                    "custom_field_config": {
                        "project_data_fields": [{
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
                        }],
                        "name": "A name set by the front end"
                    },
                    "name": "A name set by the front end"
                }



@given("I create a project JSON by adding one of these project types and some custom fields and Bar as a name")
def build_project_json(context):
    non_saved_search_project_type = ""
    for ptype in context.project_types:
        if not ptype["saved_search_project_type"]:
            non_saved_search_project_type = ptype
    context.project_json = {
        "project_type": non_saved_search_project_type,
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
