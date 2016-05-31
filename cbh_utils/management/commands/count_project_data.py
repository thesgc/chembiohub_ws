from django.core.management.base import BaseCommand, CommandError

import elasticsearch
class Command(BaseCommand):

    def handle(self, *args, **options):
        from cbh_chembl_model_extension.models import CBHCompoundBatch
        from cbh_core_model.models import Project
        from cbh_utils.elasticsearch_client import get_project_index_name

        ps = Project.objects.all()
        es = elasticsearch.Elasticsearch()
        for p in ps:
            ind = get_project_index_name(p.id)

            data = es.search(ind, body={"query" : {"match_all" :{}}}, ignore_unavailable=True)
            print "%s\t%d\t%d\t%d" % (p.name, p.id, data["hits"]["total"], CBHCompoundBatch.objects.filter(project=p).count())

