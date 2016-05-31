from django.core.management.base import BaseCommand, CommandError
from django.core.paginator import Paginator
try:
    # django >= 1.7
    from django.apps import apps
    get_model = apps.get_model
except ImportError:
    # django < 1.7
    from django.db.models import get_model
import gc

def queryset_iterator(queryset, chunksize=1000):
    '''''
    Iterate over a Django Queryset ordered by the primary key

    This method loads a maximum of chunksize (default: 1000) rows in it's
    memory at the same time while django normally would load all rows in it's
    memory. Using the iterator() method only causes it to not preload all the
    classes.

    Note that the implementation of the iterator does not support ordered query sets.
    '''
    pk = 0
    last_pk = queryset.order_by('-pk')[0].pk
    queryset = queryset.order_by('pk')

    while pk < last_pk:
        for row in queryset.filter(pk__gt=pk)[:chunksize]:
            pk = row.pk
            yield row
        gc.collect()


def add_ids_to_compounds():
    """Assign project IDs to compounds by saving them"""
    CBHCompoundBatch = get_model("cbh_chembl_model_extension", "CBHCompoundBatch")
    Project = get_model("cbh_core_model", "Project")
    for p in list(Project.objects.all()):
        cs = CBHCompoundBatch.objects.filter(project=p).exclude(project_counter=-1).order_by("id")

        if cs.count() > 0:
            for obj in queryset_iterator(cs):

                obj.save()

        print "done %s" % p.name




class Command(BaseCommand):

    def handle(self, *args, **options):
        add_ids_to_compounds()







