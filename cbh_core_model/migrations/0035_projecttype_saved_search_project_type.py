# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0034_auto_20160122_0057'),
    ]

    operations = [
        migrations.AddField(
            model_name='projecttype',
            name='saved_search_project_type',
            field=models.BooleanField(default=False),
            preserve_default=True,
        ),
    ]
