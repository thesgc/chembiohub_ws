# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0024_auto_20151203_1058'),
    ]

    operations = [
        migrations.AlterField(
            model_name='invitation',
            name='email',
            field=models.CharField(max_length=100),
            preserve_default=True,
        ),
    ]
