from django.core.management.base import BaseCommand, CommandError
from django.http import HttpRequest

class Command(BaseCommand):

    def handle(self, *args, **options):
        from cbh_chem_api.resources import IndexingCBHCompoundBatchResource

        cbr = IndexingCBHCompoundBatchResource()
        cbr.reindex_elasticsearch(HttpRequest())
