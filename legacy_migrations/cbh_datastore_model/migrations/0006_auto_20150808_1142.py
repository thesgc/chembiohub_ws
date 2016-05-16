# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model',
         '0005_datapointclassification_data_form_config'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='datapointclassification',
            unique_together=set(
                [('data_form_config', 'l0', 'l1', 'l2', 'l3', 'l4')]),
        ),
    ]
