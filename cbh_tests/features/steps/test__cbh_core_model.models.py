from behave import given, when, then
import json
from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType, PERMISSION_CODENAME_SEPARATOR
from cbh_datastore_model.models import DataPoint, DataPointClassification, DataPointClassificationPermission
from django.db import IntegrityError
from django.contrib.auth.models import User, Permission
from test__cbh_core_ws__resources import logintestuser, logout


@then("the project permission name matches the new name of the project")
def test_sync_permissions_name_change(context):

    projdata = json.loads(context.updated_project_response.content)
    id = projdata["id"]
    perms = Permission.objects.filter(codename__startswith=str(id) + PERMISSION_CODENAME_SEPARATOR)
    context.test_case.assertEqual(perms.count(), 3)
    for perm in perms:
        context.test_case.assertTrue(projdata["name"] in perm.name)


@then("the custom field config name matches the new name of the project")
def test_custom_field_config_permissions_name_change(context):
    projdata = json.loads(context.updated_project_response.content)
    context.test_case.assertTrue(projdata["name"] in projdata["custom_field_config"]["name"])

@given("I make testuser a viewer of the first project in the list")
def test_make_viewer(context):
    p = Project.objects.get(pk=context.projects_on_system[0]["id"])
    context.test_case.assertFalse(context.user.has_perm("cbh_core_model.%d%sviewer" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )))    
    p.make_viewer(context.user)
    #Need to re fetch the user from the database - see  https://docs.djangoproject.com/en/1.8/topics/auth/default/#permission-caching
    context.user = User.objects.get(pk=context.user.pk)
    context.test_case.assertTrue(context.user.has_perm("cbh_core_model.%d%sviewer" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )))


@given("I make testuser an editor of the first project in the list")
def test_make_editor(context):
    p = Project.objects.get(pk=context.projects_on_system[0]["id"])
    context.test_case.assertFalse(context.user.has_perm("cbh_core_model.%d%seditor" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )))
    context.user = p.make_editor(context.user)
    #Need to re fetch the user from the database - see  https://docs.djangoproject.com/en/1.8/topics/auth/default/#permission-caching
    context.user = User.objects.get(pk=context.user.pk)
    context.test_case.assertTrue(context.user.has_perm("cbh_core_model.%d%seditor" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )))

@given("I remove all of the testusers permissions")
def remove_perms(context):
    context.user.user_permissions.clear()
    #Need to re fetch the user from the database - see  https://docs.djangoproject.com/en/1.8/topics/auth/default/#permission-caching
    context.user = User.objects.get(pk=context.user.pk)

