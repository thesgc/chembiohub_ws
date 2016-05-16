# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model',
         '0008_datapointclassificationpermission_data_point_classification'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='datapointclassificationpermission',
            name='created_by',
        ),
        migrations.RemoveField(
            model_name='datapointclassificationpermission',
            name='data_point_classification',
        ),
        migrations.RemoveField(
            model_name='datapointclassificationpermission',
            name='project',
        ),
        migrations.DeleteModel(
            name='DataPointClassificationPermission',
        ),
    ]
