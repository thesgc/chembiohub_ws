# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0004_auto_20150327_1205'),
    ]

    operations = [
        migrations.AddField(
            model_name='pinnedcustomfield',
            name='field_key',
            field=models.CharField(default=b'', max_length=50),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='field_type',
            field=models.CharField(default=b'char', max_length=15, choices=[(b'textarea', b'Full text'), (b'text', b'Short text field'), (b'uiselecttag', b'Choice allowing create'), (
                b'number', b'Decimal field'), (b'uiselect', b'Choice field'), (b'integer', b'Integer field'), (b'uiselecttags', b'Tags field allowing create')]),
            preserve_default=True,
        ),
    ]
