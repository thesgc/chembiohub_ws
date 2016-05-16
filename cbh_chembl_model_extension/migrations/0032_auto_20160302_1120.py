# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0031_auto_20151203_0856'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='cbhcompoundbatch',
            name='editable_by',
        ),
        migrations.RemoveField(
            model_name='cbhcompoundbatch',
            name='errors',
        ),
    ]
