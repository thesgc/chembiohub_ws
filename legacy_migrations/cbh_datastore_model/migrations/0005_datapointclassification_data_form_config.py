# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0008_auto_20150808_1140'),
        ('cbh_datastore_model', '0004_auto_20150808_1139'),
    ]

    operations = [
        migrations.AddField(
            model_name='datapointclassification',
            name='data_form_config',
            field=models.ForeignKey(
                default=None, to='cbh_core_model.DataFormConfig'),
            preserve_default=False,
        ),
    ]
