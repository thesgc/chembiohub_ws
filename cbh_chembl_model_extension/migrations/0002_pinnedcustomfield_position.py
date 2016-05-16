# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='pinnedcustomfield',
            name='position',
            field=models.PositiveSmallIntegerField(default=1),
            preserve_default=False,
        ),
    ]
