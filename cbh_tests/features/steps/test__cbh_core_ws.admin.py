"""
Note I do not test the admin UI directly as this is too much work for little reward, however I do test the custom function which is called via the admin ui
"""


from behave import given, when, then
import json
from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType, PERMISSION_CODENAME_SEPARATOR
from django.db import IntegrityError
from django.contrib.auth.models import User, Permission




@when("I save the data form config via the admin ui")
def data_form_config_admin(context):
    context.dfc.get_all_ancestor_objects(context)


@then("All ancestor data form configs are present")
def ancestor_data_form_configs(context):
    context.test_case.assertEqual(context.dfc.parent.parent.parent.l0_id, context.dfc.l0_id)
    context.test_case.assertEqual(context.dfc.parent.parent.parent.l1_id, None)
    context.test_case.assertEqual(context.dfc.parent.parent.parent.l2_id, None)
    context.test_case.assertEqual(context.dfc.parent.parent.parent.l3_id, None)

    context.test_case.assertEqual(context.dfc.parent.parent.l0_id, context.dfc.l0_id)
    context.test_case.assertEqual(context.dfc.parent.parent.l1_id, context.dfc.l1_id)
    context.test_case.assertEqual(context.dfc.parent.parent.l2_id, None)
    context.test_case.assertEqual(context.dfc.parent.parent.l3_id, None)


    context.test_case.assertEqual(context.dfc.parent.l0_id, context.dfc.l0_id)
    context.test_case.assertEqual(context.dfc.parent.l1_id, context.dfc.l1_id)
    context.test_case.assertEqual(context.dfc.parent.l2_id, context.dfc.l2_id)
    context.test_case.assertEqual(context.dfc.parent.l3_id, None)




@given("I take the first project in the list and I link it to my data form config via admin ui")
def data_form_config_admin(context):
    p = Project.objects.get(pk=context.projects_on_system[0]["id"])
    p.enabled_forms.add(context.dfc)
    p.save()
    
