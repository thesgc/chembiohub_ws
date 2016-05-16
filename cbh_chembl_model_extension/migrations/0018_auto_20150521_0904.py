# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0017_auto_20150512_0831'),
    ]

    operations = [
        migrations.AddField(
            model_name='cbhcompoundbatch',
            name='batch_number',
            field=models.IntegerField(default=-1, null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cbhcompoundbatch',
            name='blinded_batch_id',
            field=models.CharField(
                default=b'', max_length=12, null=True, blank=True),
            preserve_default=True,
        ),
    ]
