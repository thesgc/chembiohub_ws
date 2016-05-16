# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


def create_data_types(apps, schema_editor):
    #This has been replaced with a fixture
    pass




class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0017_auto_20150915_0022'),
    ]

    operations = [
        migrations.RunPython(create_data_types),
    ]
