# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0007_auto_20150330_1052'),
    ]

    operations = [
        migrations.AddField(
            model_name='customfieldconfig',
            name='schemaform',
            field=models.TextField(default=b''),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='field_key',
            field=models.CharField(default=b'', unique=True, max_length=50),
            preserve_default=True,
        ),
    ]
