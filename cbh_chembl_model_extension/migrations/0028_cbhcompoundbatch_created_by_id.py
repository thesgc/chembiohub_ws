# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0027_auto_20151103_0801'),
    ]

    operations = [
        migrations.AddField(
            model_name='cbhcompoundbatch',
            name='created_by_id',
            field=models.IntegerField(default=None, null=True, blank=True),
            preserve_default=True,
        ),
    ]
