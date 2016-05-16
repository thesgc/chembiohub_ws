# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0012_auto_20150911_0825'),
    ]

    operations = [
        migrations.AddField(
            model_name='dataformconfig',
            name='human_added',
            field=models.NullBooleanField(default=True),
            preserve_default=True,
        ),
    ]
