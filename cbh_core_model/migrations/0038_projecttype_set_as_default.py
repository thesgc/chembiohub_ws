# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0037_auto_20160127_0700'),
    ]

    operations = [
        migrations.AddField(
            model_name='projecttype',
            name='set_as_default',
            field=models.BooleanField(default=False),
            preserve_default=True,
        ),
    ]
