# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model', '0014_auto_20150819_1503'),
    ]

    operations = [
        migrations.AddField(
            model_name='datapointclassification',
            name='parent',
            field=models.ForeignKey(related_name='children', default=None,
                                    to='cbh_datastore_model.DataPointClassification', null=True),
            preserve_default=True,
        ),
    ]
