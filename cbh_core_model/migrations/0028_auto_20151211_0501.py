# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0027_auto_20151211_0451'),
    ]

    operations = [
        migrations.AddField(
            model_name='skinningconfig',
            name='data_manager_email',
            field=models.CharField(default=b'', max_length=100),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='skinningconfig',
            name='data_manager_name',
            field=models.CharField(default=b'', max_length=100),
            preserve_default=True,
        ),
    ]
