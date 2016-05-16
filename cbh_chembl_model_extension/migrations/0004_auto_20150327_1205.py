# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import django_extensions.db.fields


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0003_auto_20150327_1145'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='description',
            field=models.CharField(
                default=b'', max_length=1024, null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='part_of_blinded_key',
            field=models.BooleanField(
                default=False, verbose_name=b'blind key'),
            preserve_default=True,
        ),
    ]
