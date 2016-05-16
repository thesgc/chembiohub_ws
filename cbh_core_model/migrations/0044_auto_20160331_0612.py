# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0043_auto_20160303_0537'),
    ]

    operations = [
        migrations.AlterField(
            model_name='projecttype',
            name='custom_field_config_template_id',
            field=models.IntegerField(default=None, null=True),
            preserve_default=True,
        ),
    ]
