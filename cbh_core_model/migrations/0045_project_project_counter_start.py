# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0044_auto_20160331_0612'),
    ]

    operations = [
        migrations.AddField(
            model_name='project',
            name='project_counter_start',
            field=models.IntegerField(default=1),
            preserve_default=True,
        ),
    ]
