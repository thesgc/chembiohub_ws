# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0003_pinnedcustomfield_default'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='default',
            field=models.CharField(default=b'', max_length=500, blank=True),
            preserve_default=True,
        ),
    ]
