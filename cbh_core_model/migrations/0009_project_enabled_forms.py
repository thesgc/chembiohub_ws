# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0008_auto_20150808_1140'),
    ]

    operations = [
        migrations.AddField(
            model_name='project',
            name='enabled_forms',
            field=models.ManyToManyField(to='cbh_core_model.DataFormConfig'),
            preserve_default=True,
        ),
    ]
