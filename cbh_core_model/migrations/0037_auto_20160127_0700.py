# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0036_projecttype_custom_field_config_template'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='projecttype',
            name='custom_field_config_template',
        ),
        migrations.AddField(
            model_name='projecttype',
            name='custom_field_config_template_id',
            field=models.IntegerField(default=-1),
            preserve_default=True,
        ),
    ]
