# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension',
         '0020_remove_skinningconfig_created_by'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='skinningconfig',
            options={'verbose_name': 'Skinning Configuration'},
        ),
        migrations.RemoveField(
            model_name='skinningconfig',
            name='created',
        ),
        migrations.RemoveField(
            model_name='skinningconfig',
            name='modified',
        ),
    ]
