# -*- coding: utf-8 -*-
"""environment -- environmental setup for Django+Behave+Mechanize

This should go in features/environment.py
(http://packages.python.org/behave/tutorial.html#environmental-controls)

Requirements:
http://pypi.python.org/pypi/behave/
http://pypi.python.org/pypi/mechanize/
http://pypi.python.org/pypi/wsgi_intercept
http://pypi.python.org/pypi/BeautifulSoup/

Acknowledgements:
For the basic solution: https://github.com/nathforge/django-mechanize/

"""


import os
import django
import urlparse
# This is necessary for all installed apps to be recognized, for some reason.
# Already set this on vagrant
#os.environ['DJANGO_SETTINGS_MODULE'] = 'myproject.settings'
from tastypie.test import TestApiClient
from tastypie.serializers import Serializer


def before_all(context):
    # Even though DJANGO_SETTINGS_MODULE is set, this may still be
    # necessary. Or it may be simple CYA insurance.

    # We'll use thise later to frog-march Django through the motions
    # of setting up and tearing down the test environment, including
    # test databases.
    # from django.core.management import setup_environ
    # from deployment.settings  import development as settings
    # setup_environ(settings)
    #os.environ.setdefault("DJANGO_SETTINGS_MODULE", "myapp.settings")
    django.setup()

    # Take a TestRunner hostage.
    from django.test.simple import DjangoTestSuiteRunner
    context.runner = DjangoTestSuiteRunner(interactive=False)

    # If you use South for migrations, uncomment this to monkeypatch
    # syncdb to get migrations to run.
    #from south.management.commands import patch_for_test_db_setup
    # patch_for_test_db_setup()

    context.runner.setup_test_environment()
    context.runner.setup_databases()


def before_scenario(context, scenario):
    # Set up the scenario test environment

    # Set up the WSGI intercept "port".
    context.api_client = TestApiClient()
    context.ser = Serializer()
    context.post_data = {}
    context.batch = None
    context.g = None
    context.u = None
    context.response = None
    context.runner.setup_test_environment()
    # We must set up and tear down the entire database between
    # scenarios. We can't just use db transactions, as Django's
    # TestClient does, if we're doing full-stack tests with Mechanize,
    # because Django closes the db connection after finishing the HTTP
    # response.
    # Set up the Mechanize browser.
    # from wsgi_intercept import mechanize_intercept
    # # MAGIC: All requests made by this monkeypatched browser to the magic
    # # host and port will be intercepted by wsgi_intercept via a
    # # fake socket and routed to Django's WSGI interface.
    # browser = context.browser = mechanize_intercept.Browser()
    # browser.set_handle_robots(False)


def after_scenario(context, scenario):
    # Tear down the scenario test environment.
    # context.runner.teardown_databases(context.old_db_config)
    context.api_client.client.logout()

    from cbh_chembl_model_extension.models import CBHCompoundBatch
    from cbh_core_model.models import Project, CustomFieldConfig

    from django.contrib.auth.models import User, Group
    User.objects.all().delete()
    Project.objects.all().delete()
    CustomFieldConfig.objects.all().delete()
    Group.objects.all().delete()
    CBHCompoundBatch.objects.all().delete()

    context.runner.teardown_test_environment()
    # Bob's your uncle.
