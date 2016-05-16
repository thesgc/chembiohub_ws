# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_id_generator', '0006_cbhplugin_name'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='cbhplugin',
            name='output_json_path',
        ),
    ]
