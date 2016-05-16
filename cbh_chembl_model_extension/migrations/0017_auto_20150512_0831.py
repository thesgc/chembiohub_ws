# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0016_auto_20150420_0443'),
    ]

    operations = [
        migrations.AddField(
            model_name='project',
            name='is_default',
            field=models.BooleanField(default=False),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='customfieldconfig',
            name='schemaform',
            field=models.TextField(default=b'', null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='field_type',
            field=models.CharField(default=b'char', max_length=15, choices=[(b'char', b'Short text field'), (b'uiselect', b'Choice field'), (b'textarea', b'Full text'), (b'integer', b'Integer field'), (
                b'date', b'Date Field - today or past'), (b'text', b'Short text field'), (b'percentage', b'Percentage field'), (b'uiselecttags', b'Tags field allowing create'), (b'number', b'Decimal field'), (b'uiselecttag', b'Choice allowing create')]),
            preserve_default=True,
        ),
    ]
