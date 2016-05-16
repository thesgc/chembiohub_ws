# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model', '0006_auto_20150808_1142'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='datapointclassification',
            name='l0_permission',
        ),
    ]
