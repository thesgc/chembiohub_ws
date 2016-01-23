from behave import given, when, then
import json
from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType, PERMISSION_CODENAME_SEPARATOR
from cbh_datastore_model.models import DataPoint, DataPointClassification, DataPointClassificationPermission
from django.db import IntegrityError
from django.contrib.auth.models import User, Permission

@when(u'I get my single project via the data form config API as in assayreg data overview')
def test_cbh_projects_with_forms_get(context):
    from django.conf import settings
    resp = context.api_client.client.get("/" + settings.WEBSERVICES_NAME + "/datastore/cbh_projects_with_forms/?project_key=" + context.projects_on_system[0]["project_key"])
    context.project_with_forms = resp
    context.test_case.assertHttpOK(resp)


@then(u'there is a nest of data form configs down to l3')
def step_impl(context):
    data = json.loads(context.project_with_forms.content)["objects"][0]
    print (data)
    context.test_case.assertEquals(
        data["data_form_configs"][0]["last_level"], "l0")
    context.test_case.assertEquals(len(data["data_form_configs"][0][
                                   "permitted_children"]), 1)



@given(u'I also add 2 data form configs to that project')
def step_impl(context):
    raise NotImplementedError(u'STEP: Given I also add 2 data form configs to that project')

@when(u'I POST a project to cbh_projects_with_forms')
def step_impl(context):
    raise NotImplementedError(u'STEP: When I POST a project to cbh_projects_with_forms')


