# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0026_auto_20151208_0626'),
    ]

    operations = [
        migrations.AddField(
            model_name='skinningconfig',
            name='enable_smiles_input',
            field=models.NullBooleanField(default=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='skinningconfig',
            name='file_errors_from_backend',
            field=models.NullBooleanField(default=False),
            preserve_default=True,
        ),
    ]
