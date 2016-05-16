# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model', '0002_auto_20150808_1120'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='dataformconfig',
            unique_together=set([('l0', 'l1', 'l2', 'l3', 'l4')]),
        ),
    ]
