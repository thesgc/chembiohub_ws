# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import time
import os

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
from tastypie.test import ResourceTestCase
from django.conf import settings
from django.test import TestCase
from django.test.client import FakePayload, Client

from tastypie.serializers import Serializer

try:
    from urllib.parse import urlparse
except ImportError:
    from urlparse import urlparse
import unittest


class Tester(unittest.TestCase):
    # def assertEqual(a,b):
    #     return assert(a==b)

    def runTest(self):
        pass

    def assertHttpOK(self, resp):
        """
        Ensures the response is returning a HTTP 200.
        """
        return self.assertEqual(resp.status_code, 200)

    def assertHttpCreated(self, resp):
        """
        Ensures the response is returning a HTTP 201.
        """
        return self.assertEqual(resp.status_code, 201)

    def assertHttpAccepted(self, resp):
        """
        Ensures the response is returning either a HTTP 202 or a HTTP 204.
        """
        return self.assertIn(resp.status_code, [202, 204])

    def assertHttpMultipleChoices(self, resp):
        """
        Ensures the response is returning a HTTP 300.
        """
        return self.assertEqual(resp.status_code, 300)

    def assertHttpSeeOther(self, resp):
        """
        Ensures the response is returning a HTTP 303.
        """
        return self.assertEqual(resp.status_code, 303)

    def assertHttpNotModified(self, resp):
        """
        Ensures the response is returning a HTTP 304.
        """
        return self.assertEqual(resp.status_code, 304)

    def assertHttpBadRequest(self, resp):
        """
        Ensures the response is returning a HTTP 400.
        """
        return self.assertEqual(resp.status_code, 400)

    def assertHttpUnauthorized(self, resp):
        """
        Ensures the response is returning a HTTP 401.
        """
        return self.assertEqual(resp.status_code, 401)

    def assertHttpForbidden(self, resp):
        """
        Ensures the response is returning a HTTP 403.
        """
        return self.assertEqual(resp.status_code, 403)

    def assertHttpNotFound(self, resp):
        """
        Ensures the response is returning a HTTP 404.
        """
        return self.assertEqual(resp.status_code, 404)

    def assertHttpMethodNotAllowed(self, resp):
        """
        Ensures the response is returning a HTTP 405.
        """
        return self.assertEqual(resp.status_code, 405)

    def assertHttpConflict(self, resp):
        """
        Ensures the response is returning a HTTP 409.
        """
        return self.assertEqual(resp.status_code, 409)

    def assertHttpGone(self, resp):
        """
        Ensures the response is returning a HTTP 410.
        """
        return self.assertEqual(resp.status_code, 410)

    def assertHttpUnprocessableEntity(self, resp):
        """
        Ensures the response is returning a HTTP 422.
        """
        return self.assertEqual(resp.status_code, 422)

    def assertHttpTooManyRequests(self, resp):
        """
        Ensures the response is returning a HTTP 429.
        """
        return self.assertEqual(resp.status_code, 429)

    def assertHttpApplicationError(self, resp):
        """
        Ensures the response is returning a HTTP 500.
        """
        return self.assertEqual(resp.status_code, 500)

    def assertHttpNotImplemented(self, resp):
        """
        Ensures the response is returning a HTTP 501.
        """
        return self.assertEqual(resp.status_code, 501)

    def assertValidJSON(self, data):
        """
        Given the provided ``data`` as a string, ensures that it is valid JSON &
        can be loaded properly.
        """
        # Just try the load. If it throws an exception, the test case will
        # fail.
        import json
        json.loads(data)

    def assertValidXML(self, data):
        """
        Given the provided ``data`` as a string, ensures that it is valid XML &
        can be loaded properly.
        """
        # Just try the load. If it throws an exception, the test case will
        # fail.
        self.serializer.from_xml(data)

    def assertValidYAML(self, data):
        """
        Given the provided ``data`` as a string, ensures that it is valid YAML &
        can be loaded properly.
        """
        # Just try the load. If it throws an exception, the test case will
        # fail.
        self.serializer.from_yaml(data)

    def assertValidPlist(self, data):
        """
        Given the provided ``data`` as a string, ensures that it is valid
        binary plist & can be loaded properly.
        """
        # Just try the load. If it throws an exception, the test case will
        # fail.
        self.serializer.from_plist(data)

    def assertValidJSONResponse(self, resp):
        """
        Given a ``HttpResponse`` coming back from using the ``client``, assert that
        you get back:
        * An HTTP 200
        * The correct content-type (``application/json``)
        * The content is valid JSON
        """
        self.assertHttpOK(resp)
        self.assertTrue(resp['Content-Type'].startswith('application/json'))
        from django.utils.encoding import force_text

        self.assertValidJSON(force_text(resp.content))

    def assertValidXMLResponse(self, resp):
        """
        Given a ``HttpResponse`` coming back from using the ``client``, assert that
        you get back:
        * An HTTP 200
        * The correct content-type (``application/xml``)
        * The content is valid XML
        """
        self.assertHttpOK(resp)
        self.assertTrue(resp['Content-Type'].startswith('application/xml'))
        from django.utils.encoding import force_text

        self.assertValidXML(force_text(resp.content))

    def assertValidYAMLResponse(self, resp):
        """
        Given a ``HttpResponse`` coming back from using the ``client``, assert that
        you get back:
        * An HTTP 200
        * The correct content-type (``text/yaml``)
        * The content is valid YAML
        """
        self.assertHttpOK(resp)
        self.assertTrue(resp['Content-Type'].startswith('text/yaml'))
        from django.utils.encoding import force_text

        self.assertValidYAML(force_text(resp.content))

    def assertValidPlistResponse(self, resp):
        """
        Given a ``HttpResponse`` coming back from using the ``client``, assert that
        you get back:
        * An HTTP 200
        * The correct content-type (``application/x-plist``)
        * The content is valid binary plist data
        """
        self.assertHttpOK(resp)
        self.assertTrue(resp['Content-Type'].startswith('application/x-plist'))
        from django.utils.encoding import force_text

        self.assertValidPlist(force_text(resp.content))

import os
CONDA_ENV_PATH = os.getenv("CONDA_ENV_PATH")
host = CONDA_ENV_PATH + "/var/postgressocket"
database_path = CONDA_ENV_PATH + "/var/"


class TClient(TestApiClient):

    def get_content_type(self, short_format):
        """
        Given a short name (such as ``json`` or ``xml``), returns the full
        content-type for it (``application/json`` or ``application/xml`` in
        this case).
        """
        print (short_format)
        return self.serializer.content_types.get(short_format, 'json')

def before_all(context):

    from subprocess import Popen, PIPE, call
    process = Popen(['python', 'manage.py', 'migrate', '--list'], stdout=PIPE)
    commit_based_filename, error = process.communicate()

    import hashlib
    h = hashlib.md5(commit_based_filename.encode())
    commit_based_filename = h.hexdigest()
    
    commit_based_filename = database_path + commit_based_filename
    context.commit_based_filename = commit_based_filename
    print(commit_based_filename)
    context.has_run_some_scenarios = False
    if not os.path.isfile(commit_based_filename):
        call(
            'psql -h %s -c "create extension if not exists hstore;create extension if not exists rdkit;" template1' % host , shell=True)
        call(
            "dropdb dev_db --if-exists -h %s" % host , shell=True)
        call(
            "createdb dev_db -h %s -T template1" % host, shell=True)
        call(
            "python manage.py migrate", shell=True)   
        call(
            "pg_dump dev_db -Fc -h %s > %s" % (host, commit_based_filename) , shell=True)    
            
        pass
    # Even though DJANGO_SETTINGS_MODULE is set, this may still be
    # necessary. Or it may be simple CYA insurance.

    # We'll use thise later to frog-march Django through the motions
    # of setting up and tearing down the test environment, including
    # test databases.
    # from django.core.management import setup_environ
    # from deployment.settings  import development as settings
    # setup_environ(settings)
    #os.environ.setdefault("DJANGO_SETTINGS_MODULE", "myapp.settings")

    # context.runner.setup_databases()
    # from django.core.management import call_command
    # # from django.contrib.contenttypes.models import ContentType
    # # ContentType.objects.all().delete()
    # from cbh_core_model.models import Project, CustomFieldConfig
    # # Project.objects.all().delete()
    # call_command("loaddata", "/home/vagrant/chembiohub_ws/src/cbh_datastore_ws/cbh_datastore_ws/features/fixtures/newtestfixture.json")

import time
def before_scenario(context, scenario):
    # Set up the scenario test environment
    from subprocess import Popen, PIPE, call
    context.has_run_some_scenarios = True


    call(
            "dropdb dev_db --if-exists -h %s" % host , shell=True)
    call(
            "createdb dev_db -h %s -T template1" % host, shell=True)

    call(
        "pg_restore -Fc -h %s -d dev_db < %s" % (host, context.commit_based_filename), shell=True)

    django.setup()


    from django.core.management import call_command
    #must index the data


    try:
        from django.test.simple import DjangoTestSuiteRunner
    except ImportError:
        from django.test.runner import DiscoverRunner as DjangoTestSuiteRunner
    context.runner = DjangoTestSuiteRunner(interactive=False)

    context.api_client = TClient()
    context.test_case = Tester()

    from cbh_chembl_model_extension.models import CBHCompoundBatch
    from cbh_core_model.models import Project, CustomFieldConfig, SkinningConfig
    from django.contrib.auth.models import User, Group
    context.dfc =None
    context.response = None
    context.user ,created  = User.objects.get_or_create(username="testuser")
    context.user.set_password("testuser")
    context.user.save()
    context.superuser ,created = User.objects.get_or_create(username="superuser")
    context.superuser.set_password("superuser")
    context.superuser.is_superuser = True
    context.superuser.is_staff = True
    context.superuser.save()
    from django.test import Client
    context.dclient = Client()
    context.dclient.login(username="testuser", password="testuser")
    call_command("reindex_compounds_new")   
    import elasticsearch
    es = elasticsearch.Elasticsearch()
    #Delete all existing indices
    es.indices.delete("dev*")
    context.runner.setup_test_environment()
    

def after_scenario(context, scenario):
    # Tear down the scenario test environment.
    # context.runner.teardown_databases(context.old_db_config)
    try:
        os.remove("qclusterprocesslog.log")
    except OSError:
        pass
    context.api_client.client.logout()
    
    from django.db import connection
    connection.close()



def after_all(context):
    from django.db import connection
    connection.close()
    from subprocess import Popen, PIPE, call

        #We can't import resources unless django.setup has been called
    before_scenario(context, None)
    
    from django.core.management import call_command
    call_command("reindex_compounds_new")
    
    call_command("loaddata", "projecttypes.json")
    # call_command("loaddata", "project_setup.json")


