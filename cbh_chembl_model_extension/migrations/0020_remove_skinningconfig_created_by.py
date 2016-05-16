# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0019_auto_20150721_0515'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='skinningconfig',
            name='created_by',
        ),
    ]
