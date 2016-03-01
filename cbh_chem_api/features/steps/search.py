from behave import given, when, then
import json
import os
from django.contrib.auth.models import User, Group


@when("I specify a {start_or_end} date")
def step(context, start_or_end, date):
    if start_or_end == "start":
        context.post_data[created__gte] = date(2015, 1, 6).ctime()
    elif start_or_end == "end":
        context.post_data[created__lte] = date(2015, 6, 1).ctime()


@when("I specify a project {project}")
def step(context, project):
    context.post_data[project__project_key] = project


@when("I specify a custom field tag {custom_tag}")
def step(context, custom_tag):
    context.post_data.custom_fields[custom_tag] = custom_tag


@when("I specify structural search type {struc_search_type}")
def step(context, struc_search_type):
    context.post_data[struc_search_type] = context.post_data["ctab"]


@then("I submit my search")
def step(context):
    path = "/devapi/cbh_compound_batches/"
    func = context.api_client.post
    #context.post_data["projectKey"] = projkey
    # print(context.post_data)
    resp = func(
        path,
        format='json',
        data=context.post_data,
    )
    context.returned_data = resp
    assert resp.status_code == 201


@then("I can see search results")
def step(context):
    found_cmpds = context.ser.deserialize(
        context.returned_data.content)["objects"]
    # retrieve registered inchi
    #reg_inchi = found_cmpds[0]['standardInchi']
    assert len(found_cmpds) > 0
