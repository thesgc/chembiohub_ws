from django.core.management.base import BaseCommand, CommandError
from django.contrib.sessions.models import Session


class Command(BaseCommand):

    def handle(self, *args, **options):
        """Log the users out when there has been a change to the code"""
        Session.objects.all().delete()


