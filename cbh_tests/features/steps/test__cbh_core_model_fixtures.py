from behave import given, when, then
import json
from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType, PERMISSION_CODENAME_SEPARATOR
from django.db import IntegrityError
from django.contrib.auth.models import User, Permission
from django.core.management import call_command

@given("I have loaded the fixtures for project types and data types")
def load_fixtures(context):
    call_command("loaddata", "datatypes.json")
    call_command("loaddata", "projecttypes.json")