# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model', '0011_auto_20150810_1258'),
    ]

    operations = [
        migrations.AlterField(
            model_name='datapointclassification',
            name='l0',
            field=models.ForeignKey(
                related_name='l0', default=1, to='cbh_datastore_model.DataPoint'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='datapointclassification',
            name='l1',
            field=models.ForeignKey(
                related_name='l1', default=1, to='cbh_datastore_model.DataPoint'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='datapointclassification',
            name='l2',
            field=models.ForeignKey(
                related_name='l2', default=1, to='cbh_datastore_model.DataPoint'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='datapointclassification',
            name='l3',
            field=models.ForeignKey(
                related_name='l3', default=1, to='cbh_datastore_model.DataPoint'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='datapointclassification',
            name='l4',
            field=models.ForeignKey(
                related_name='l4', default=1, to='cbh_datastore_model.DataPoint'),
            preserve_default=True,
        ),
    ]
