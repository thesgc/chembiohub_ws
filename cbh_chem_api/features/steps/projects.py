from behave import given, when, then
import json
from django.contrib.auth.models import User, Group


@then("I get projects and the response code will be {responsecode} with {project} in it")
def step(context, responsecode, project):

    resp = context.api_client.get("/dev/cbh_projects/",  format='json')
    print(resp.status_code)
    assert int(resp.status_code) == int(responsecode)
    assert project in resp.content


@then("I get projects and the response code will be {responsecode} without {project} in it")
def step(context, responsecode, project):

    resp = context.api_client.get("/dev/cbh_projects/",  format='json')
    print(resp.status_code)
    assert int(resp.status_code) == int(responsecode)
    assert project not in resp.content


@then("I get projects and the response code will be {responsecode}")
def step(context, responsecode):
    resp = context.api_client.get("/dev/cbh_projects/",  format='json')
    assert int(resp.status_code) == int(responsecode)


@given("my user is member of a group")
def step(context):
    context.g = Group.objects.create(name="g")
    context.g.user_set.add(context.u)


@given("a valid project exists {projkey}")
def step(context, projkey):
    from cbh_core_model.models import Project, CustomFieldConfig
    cfc = CustomFieldConfig.objects.create(created_by=context.u, name='test')
    Project.objects.create(
        name=projkey, project_key=projkey, created_by=context.u, custom_field_config=cfc)


@given("I automatically have editor permissions as creator")
def step(context):
    from cbh_core_model.models import Project
    for p in Project.objects.all():
        assert p.created_by.has_perm("%d.%s" % (p.id, "editor"))


@given("I remove my permissions")
def step(context):
    for perm in context.u.user_permissions.all():
        context.u.user_permissions.remove(perm)
        context.u = User.objects.get(pk=context.u.id)
    from cbh_core_model.models import Project
    for p in Project.objects.all():
        assert p.created_by.has_perm("%d.%s" % (p.id, "editor")) == False


@given("I have given {me_or_group} editor privileges for {projkey}")
def step(context, me_or_group, projkey):
    from cbh_core_model.models import Project
    p = Project.objects.get(project_key=projkey)
    if me_or_group == "myself":
        p.make_editor(context.u)
    elif me_or_group == "mygroup":
        p.make_editor(context.g)
    context.u = User.objects.get(pk=context.u.id)
    assert context.u.has_perm("%d.%s" % (p.id, "editor"))


@given("I have given {me_or_group} viewer privileges for {projkey}")
def step(context, me_or_group, projkey):
    from cbh_core_model.models import Project
    p = Project.objects.get(project_key=projkey)
    if me_or_group == "myself":
        p.make_viewer(context.u)
    elif me_or_group == "mygroup":
        p.make_viewer(context.g)
    context.u = User.objects.get(pk=context.u.id)
    assert context.u.has_perm("%d.%s" % (p.id, "viewer"))


@given("I havent given {grouporuser} {editororviewer} privileges for {projkey}")
def step(context, grouporuser, editororviewer, projkey):
    pass
