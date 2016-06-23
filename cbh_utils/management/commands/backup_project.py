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
            print ("Pick a project or comma separated list of projects to back up")
            
            for p in ps:

                print ("%d:\t%s" % ( p.id, p.name,))
            list_of_pids = input("Enter a number or list of numbers from above\n")
            pid_numbers = [int(pid.strip()) for pid in list_of_pids.split(",")]
            message = "\n\nReally backup these projects?\n"
            to_delete = []
            for pid in pid_numbers:
                for proj in ps:
                    if pid == proj.id:
                        to_delete.append(proj)
                        message += "%d:\t%s" % ( proj.id, proj.name,)
                        message += " (and %d compound batches)\n" % CBHCompoundBatch.objects.filter(project=proj).count()
                    
            message += "Type BACKUP and hit enter to say yes\n"
            delete = input(message)
            if delete == "BACKUP":
                backup_permissions(to_delete, directory)
                backup_attachments(to_delete, directory)
                backup_projects(to_delete, directory)
                backup_compounds(to_delete, directory)
                
            else:
                print ("BACKUP not confirmed, exiting")

            # for p in ps:
            #     ind = get_project_index_name(p.id)

            #     data = es.search(ind, body={"query" : {"match_all" :{}}}, ignore_unavailable=True)
            #     print "%s\t%d\t%d\t%d" % (p.name, p.id, data["hits"]["total"], CBHCompoundBatch.objects.filter(project=p).count())
        else:
            print("No projects to BACKUP, exiting.")