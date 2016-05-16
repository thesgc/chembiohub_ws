# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_id_generator', '0002_cbhcompoundid_current_batch_id'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='cbhcompoundid',
            name='original_project_key',
        ),
    ]
