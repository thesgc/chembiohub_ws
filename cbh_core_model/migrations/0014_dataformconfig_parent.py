# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0013_dataformconfig_human_added'),
    ]

    operations = [
        migrations.AddField(
            model_name='dataformconfig',
            name='parent',
            field=models.ForeignKey(
                default=None, blank=True, to='cbh_core_model.DataFormConfig', null=True),
            preserve_default=True,
        ),
    ]
