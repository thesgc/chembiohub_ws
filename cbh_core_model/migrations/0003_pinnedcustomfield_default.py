# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0002_auto_20150806_0658'),
    ]

    operations = [
        migrations.AddField(
            model_name='pinnedcustomfield',
            name='default',
            field=models.CharField(default=b'', max_length=50, blank=True),
            preserve_default=True,
        ),
    ]
