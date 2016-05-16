# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model', '0003_auto_20150808_1129'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='dataformconfig',
            unique_together=None,
        ),
        migrations.RemoveField(
            model_name='dataformconfig',
            name='created_by',
        ),
        migrations.RemoveField(
            model_name='dataformconfig',
            name='l0',
        ),
        migrations.RemoveField(
            model_name='dataformconfig',
            name='l1',
        ),
        migrations.RemoveField(
            model_name='dataformconfig',
            name='l2',
        ),
        migrations.RemoveField(
            model_name='dataformconfig',
            name='l3',
        ),
        migrations.RemoveField(
            model_name='dataformconfig',
            name='l4',
        ),
        migrations.AlterUniqueTogether(
            name='datapointclassification',
            unique_together=set([]),
        ),
        migrations.RemoveField(
            model_name='datapointclassification',
            name='data_form_config',
        ),
        migrations.DeleteModel(
            name='DataFormConfig',
        ),
    ]
