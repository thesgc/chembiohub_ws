# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0013_auto_20150409_1801'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='custom_field_config',
            field=models.ForeignKey(
                related_name='pinned_custom_field', to='cbh_chembl_model_extension.CustomFieldConfig'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='field_type',
            field=models.CharField(default=b'char', max_length=15, choices=[(b'uiselect', b'Choice field'), (b'textarea', b'Full text'), (b'integer', b'Integer field'), (b'date', b'Date Field - today or past'), (
                b'text', b'Short text field'), (b'percentage', b'Percentage field'), (b'uiselecttags', b'Tags field allowing create'), (b'number', b'Decimal field'), (b'uiselecttag', b'Choice allowing create')]),
            preserve_default=True,
        ),
    ]
