from behave import given, when, then
import json
from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType, PERMISSION_CODENAME_SEPARATOR, DataFormConfig
from cbh_datastore_model.models import DataPoint, DataPointClassification, DataPointClassificationPermission
from django.db import IntegrityError
from django.contrib.auth.models import User, Permission

@when(u'I get my single project via the data form config API as in assayreg data overview')
def test_cbh_projects_with_forms_get(context):
    from django.conf import settings
    resp = context.api_client.get("/" + settings.WEBSERVICES_NAME + "/datastore/cbh_projects_with_forms/?project_key=" + context.projects_on_system[0]["project_key"])
    context.project_with_forms = resp
    context.test_case.assertHttpOK(resp)


@given(u'I set up a project and data form config as before')
def step(context):
    context.execute_steps(u"""
        Given I have created a data form config and a project as before and I list the projects
        and I take the first project in the list and I link it to my data form config via admin ui
        When I save the data form config via the admin ui
        Then All ancestor data form configs are present
        """)


@then(u'there is a nest of data form configs down to l3')
def step_impl(context):
    data = json.loads(context.project_with_forms.content)["objects"][0]
    print (data)
    context.test_case.assertEquals(
        data["data_form_configs"][0]["last_level"], "l0")
    context.test_case.assertEquals(len(data["data_form_configs"][0][
                                   "permitted_children"]), 1)


@given("I have created a data form config and a project as before and I list the projects")
def run_steps(context):
    context.execute_steps(u"""
            Given I have loaded the fixtures for project types and data types
            Given testuser
            When I log in testuser

            Then I can list the data types on the system and there are 4
            Given I POST 4 custom field configs
            When I list the custom field configs on the system there are 5
            Given I create a data form config using the 4 custom field configs
            When I post the template data form config
            Then The data form config is created

             Given testuser has the cbh_core_model.add_project permission
            Then I can list the projects types on the system and there are 3
            Given I create a project JSON by adding one of these project types and some custom fields and Bar as a name
            When I POST a project to cbh_projects
            Then the project is created
            Then I can list the projects on the system

        """)



@when("I list the custom field configs on the system there are 5")
def step_impl(context):
    from django.conf import settings
    for postresp in context.custom_field_config_responses:
        context.test_case.assertHttpCreated(postresp)
    resp = context.api_client.get("/" + settings.WEBSERVICES_NAME + "/datastore/cbh_custom_field_config/")
    jsond = json.loads(resp.content)
    context.cfcs = jsond["objects"]
    context.test_case.assertEquals(len(context.cfcs), 5)


@given("I create a data form config using the 4 custom field configs")
def create_data_form_config(context):
    created_cfcs = [json.loads(cfc.content) for cfc in context.custom_field_config_responses]
    mapped_uris = {cfc["name"]: cfc["resource_uri"] for cfc in created_cfcs}
    mappings = (("Project", "l0"), ("Sub-Project", "l1"), ("Assay", "l2"), ("Activity", "l3"))
    
    context.data_form_config_template = {level : mapped_uris[name]  for name, level in mappings}
    print (context.data_form_config_template)

@when("I post the template data form config")
def data_form_config_post(context):
    from django.conf import settings

    resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/datastore/cbh_data_form_config/" , data=context.data_form_config_template)
    context.dfc_resp = resp


@then("The data form config is created")
def data_form_config_assert_created(context):
    context.test_case.assertHttpCreated(context.dfc_resp)
    context.dfc = DataFormConfig.objects.get(id=json.loads(context.dfc_resp.content)["id"])


@given(u'I POST 4 custom field configs')
def step_impl(context):
    from django.conf import settings

    custom_field_configs = [ {
                "project_data_fields": [{
                    "required": False,
                    "field_type": "char",
                    "open_or_restricted": "open",
                    "name": "Test Field"
                }],
                "name": "Project"
            },
            {
                "project_data_fields": [{
                    "required": False,
                    "field_type": "char",
                    "open_or_restricted": "open",
                    "name": "Test Field"
                }],
                "name": "Sub-Project"
            },
            {
                "project_data_fields": [{
                    "required": False,
                    "field_type": "char",
                    "open_or_restricted": "open",
                    "name": "Test Field"
                }],
                "name": "Assay"
            },
            {
                "project_data_fields": [{
                    "required": False,
                    "field_type": "char",
                    "open_or_restricted": "open",
                    "name": "Test Field"
                }],
                "name": "Activity"
            },]
    context.custom_field_config_responses = []
    for index, data_type in enumerate(context.data_types):
        custom_field_configs[index]["data_type"] = data_type["resource_uri"]
        resp = context.api_client.post("/" + settings.WEBSERVICES_NAME + "/datastore/cbh_custom_field_config/" , data=custom_field_configs[index])
        context.custom_field_config_responses.append(resp)





@given(u'I also add 2 data form configs to that project')
def step_impl(context):
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
            "name": "Foo",
            "data_form_configs" : [
                {

                }
            ]
        }
    raise NotImplementedError(u'STEP: Given I also add 2 data form configs to that project')

@when(u'I POST a project to cbh_projects_with_forms')
def step_impl(context):
    raise NotImplementedError(u'STEP: When I POST a project to cbh_projects_with_forms')


