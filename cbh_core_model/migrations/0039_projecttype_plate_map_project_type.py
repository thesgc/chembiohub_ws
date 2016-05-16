# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0038_projecttype_set_as_default'),
    ]

    operations = [
        migrations.AddField(
            model_name='projecttype',
            name='plate_map_project_type',
            field=models.BooleanField(default=False),
            preserve_default=True,
        ),
    ]
