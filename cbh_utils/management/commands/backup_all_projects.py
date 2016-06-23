from django.core.management.base import BaseCommand, CommandError
try: input = raw_input
except NameError: pass
import elasticsearch
import os
from django.contrib.auth import get_user_model
from django.test import RequestFactory
import shutil
from cbh_utils import elasticsearch_client
from cbh_core_api.tasks import remove_session_cached_projectlists

from cbh_utils.management.utils import  backup_projects, backup_compounds, backup_filename, backup_attachments, backup_permissions, delete_projects




class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('directory')

    def handle(self, *args, **options):
        directory = options["directory"]
        if not os.path.exists(directory):
            os.makedirs(directory)
            print "created directory"
        from cbh_chembl_model_extension.models import CBHCompoundBatch
        from cbh_core_model.models import Project
        from cbh_utils.elasticsearch_client import get_project_index_name

        ps = list(Project.objects.all().order_by("id"))
        es = elasticsearch.Elasticsearch()
        if len(ps) > 0:
            backup_permissions(ps, directory)
            backup_attachments(ps, directory)
            backup_projects(ps, directory)
            backup_compounds(ps, directory)
        else:
            print("No projects to BACKUP, exiting.")