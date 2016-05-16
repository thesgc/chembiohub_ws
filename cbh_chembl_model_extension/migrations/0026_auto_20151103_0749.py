# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0025_auto_20150807_0933'),
    ]

    operations = [
        migrations.AddField(
            model_name='cbhcompoundbatch',
            name='big_image',
            field=models.TextField(default=b''),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cbhcompoundbatch',
            name='image',
            field=models.TextField(default=b''),
            preserve_default=True,
        ),
    ]
