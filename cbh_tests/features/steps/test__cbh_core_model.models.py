from behave import given, when, then
import json

from django.db import IntegrityError


@then("the project permission name matches the new name of the project")
def test_sync_permissions_name_change(context):
    from django.contrib.auth.models import User, Permission
    from cbh_core_model.models import  PERMISSION_CODENAME_SEPARATOR
    projdata = json.loads(context.updated_project_response.content)
    id = projdata["id"]
    perms = Permission.objects.filter(codename__startswith=str(id) + PERMISSION_CODENAME_SEPARATOR)
    context.test_case.assertEqual(perms.count(), 3)
    for perm in perms:
        context.test_case.assertTrue("Foo" in perm.name)



@then(u'the project key has changed')
def step_impl(context):
    projdata = json.loads(context.updated_project_response.content)
    context.test_case.assertEqual(projdata["project_key"], "foo")




@then("the custom field config name matches the new name of the project")
def test_custom_field_config_permissions_name_change(context):
    projdata = json.loads(context.updated_project_response.content)
    context.test_case.assertTrue("Foo" in projdata["custom_field_config"]["name"])

@then("the project creator is automatically an owner")
def test_sync_permissions_owner(context):
    from cbh_core_model.models import User
    from cbh_core_model.models import  PERMISSION_CODENAME_SEPARATOR
    context.user = User.objects.get(pk=context.user.pk)
    context.test_case.assertTrue(context.user.has_perm("cbh_core_model.%d%sowner" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )))


@given("I make testuser a viewer of the first project in the list")
def test_make_viewer(context):
    from cbh_core_model.models import  PERMISSION_CODENAME_SEPARATOR
    from cbh_core_model.models import User
    from cbh_core_model.models import User, Project
    p = Project.objects.get(pk=context.projects_on_system[0]["id"])
    context.test_case.assertFalse(context.user.has_perm("cbh_core_model.%d%sviewer" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )))    
    p.make_viewer(context.user)
    #Need to re fetch the user from the database - see  https://docs.djangoproject.com/en/1.8/topics/auth/default/#permission-caching
    context.user = User.objects.get(pk=context.user.pk)
    context.test_case.assertTrue(context.user.has_perm("cbh_core_model.%d%sviewer" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )))


@given("I make testuser an editor of the first project in the list")
def test_make_editor(context):
    from cbh_core_model.models import  PERMISSION_CODENAME_SEPARATOR
    from cbh_core_model.models import User, Project
    p = Project.objects.get(pk=context.projects_on_system[0]["id"])
    context.test_case.assertFalse(context.user.has_perm("cbh_core_model.%d%seditor" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )))
    p.make_editor(context.user)
    #Need to re fetch the user from the database - see  https://docs.djangoproject.com/en/1.8/topics/auth/default/#permission-caching
    context.user = User.objects.get(pk=context.user.pk)
    context.test_case.assertTrue(context.user.has_perm("cbh_core_model.%d%seditor" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )))



@given("I make testuser an owner of the first project in the list")
def step(context):
    from cbh_core_model.models import  PERMISSION_CODENAME_SEPARATOR
    from cbh_core_model.models import User, Project
    p = Project.objects.get(pk=context.projects_on_system[0]["id"])
    #context.test_case.assertFalse(context.user.has_perm("cbh_core_model.%d%sowner" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )))
    p.make_owner(context.user)
    #Need to re fetch the user from the database - see  https://docs.djangoproject.com/en/1.8/topics/auth/default/#permission-caching
    context.user = User.objects.get(pk=context.user.pk)
    context.test_case.assertTrue(context.user.has_perm("cbh_core_model.%d%sowner" % (context.projects_on_system[0]["id"],  PERMISSION_CODENAME_SEPARATOR  )))







@given("I remove all of the testusers permissions")
def remove_perms(context):
    from cbh_core_model.models import User
    context.user.user_permissions.clear()
    #Need to re fetch the user from the database - see  https://docs.djangoproject.com/en/1.8/topics/auth/default/#permission-caching
    context.user = User.objects.get(pk=context.user.pk)






    