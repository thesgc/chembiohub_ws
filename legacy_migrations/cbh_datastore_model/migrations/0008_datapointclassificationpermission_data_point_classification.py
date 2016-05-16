# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model',
         '0007_remove_datapointclassification_l0_permission'),
    ]

    operations = [
        migrations.AddField(
            model_name='datapointclassificationpermission',
            name='data_point_classification',
            field=models.ForeignKey(related_name='l0_permission', default=None,
                                    to='cbh_datastore_model.DataPointClassificationPermission'),
            preserve_default=False,
        ),
    ]
