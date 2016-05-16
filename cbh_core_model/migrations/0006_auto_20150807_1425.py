# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


def create_data_types(apps, schema_editor):
    #This method was removed in favour of adding a fixture
    pass


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0005_auto_20150807_1425'),
    ]

    operations = [
        migrations.RunPython(create_data_types),
    ]
