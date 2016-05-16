# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0010_auto_20150809_0021'),
    ]

    operations = [
        migrations.AlterField(
            model_name='customfieldconfig',
            name='name',
            field=models.CharField(unique=True, max_length=100),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='project',
            name='name',
            field=models.CharField(
                default=None, max_length=100, null=True, db_index=True, blank=True),
            preserve_default=True,
        ),
    ]
