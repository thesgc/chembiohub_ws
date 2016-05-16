# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_id_generator', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='cbhcompoundid',
            name='current_batch_id',
            field=models.IntegerField(default=0),
            preserve_default=True,
        ),
    ]
