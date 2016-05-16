# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0018_auto_20150915_1238'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='field_key',
            field=models.CharField(default=b'', max_length=500),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='name',
            field=models.CharField(max_length=100),
            preserve_default=True,
        ),
    ]
