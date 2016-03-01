# -*- coding: utf-8 -*-
"""steps/browser_steps.py -- step implementation for our browser feature demonstration.
"""
from behave import given, when, then
import json


@given('a user')
def step(context):
    from django.contrib.auth.models import User
    try:
        u = User(username='foo', email='foo@example.com')
        u.set_password('bar')
        u.save()
    except:
        pass
    context.u = u


@when('I do not log in')
def step(context):
    pass


@then('I see a 401 error')
def step(context):

    resp = context.api_client.get("/dev/users/",  format='json')

    assert resp.status_code == 401
