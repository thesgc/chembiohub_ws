# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0016_auto_20150914_2301'),
    ]

    operations = [
        migrations.AlterField(
            model_name='datatype',
            name='name',
            field=models.CharField(unique=True, max_length=100),
            preserve_default=True,
        ),
    ]
