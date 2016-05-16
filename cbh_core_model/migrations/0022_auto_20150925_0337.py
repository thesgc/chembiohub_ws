# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model',
         '0021_pinnedcustomfield_attachment_field_mapped_to'),
    ]

    operations = [
        migrations.AlterField(
            model_name='customfieldconfig',
            name='name',
            field=models.CharField(unique=True, max_length=500),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='datatype',
            name='name',
            field=models.CharField(unique=True, max_length=500),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='name',
            field=models.CharField(max_length=500),
            preserve_default=True,
        ),
    ]
