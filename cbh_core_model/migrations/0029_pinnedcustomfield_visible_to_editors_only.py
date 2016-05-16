# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0028_auto_20151211_0501'),
    ]

    operations = [
        migrations.AddField(
            model_name='pinnedcustomfield',
            name='visible_to_editors_only',
            field=models.BooleanField(default=False),
            preserve_default=True,
        ),
    ]
