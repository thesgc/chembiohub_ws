from django.core.management.base import BaseCommand, CommandError


class Command(BaseCommand):

    def handle(self, *args, **options):
        from cbh_chembl_ws_extension.parser import ChemblAPIConverter

        ChemblAPIConverter().write_schema()
