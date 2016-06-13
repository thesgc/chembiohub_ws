# -*- coding: utf-8 -*-
# Generated by Django 1.9 on 2016-06-09 09:30
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0050_auto_20160601_0948'),
    ]

    operations = [
        migrations.AddField(
            model_name='skinningconfig',
            name='max_chem_upload_size',
            field=models.IntegerField(default=5000),
        ),
        migrations.AddField(
            model_name='skinningconfig',
            name='max_non_chem_upload_size',
            field=models.IntegerField(default=50000),
        ),
    ]
