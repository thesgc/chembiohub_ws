# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0002_pinnedcustomfield_position'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='pinnedcustomfield',
            options={'ordering': ['position'], 'get_latest_by': 'created'},
        ),
        migrations.AddField(
            model_name='pinnedcustomfield',
            name='description',
            field=models.TextField(default=b''),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='field_type',
            field=models.CharField(default=b'char', max_length=15, choices=[(b'text', b'Short text field'), (b'textarea', b'Full text'), (b'uiselect', b'Choice field'), (
                b'integer', b'Integer field'), (b'number', b'Decimal field'), (b'uiselecttag', b'Choice allowing create'), (b'uiselecttags', b'Tags field allowing create')]),
            preserve_default=True,
        ),
    ]
