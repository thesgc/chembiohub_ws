# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations

   

class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0041_auto_20160303_0441'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='field_type',
            field=models.CharField(default=b'char', max_length=15, choices=[(b'text', b'Short text field'), (b'char', b'Short text field'), (b'textarea', b'Full text'), (b'uiselect', b'Choice field'), (b'integer', b'Integer field'), (b'number', b'Decimal field'), (b'uiselecttag', b'Choice allowing create'), (b'uiselecttags', b'Tags field allowing create'), (b'percentage', b'Percentage field'), (b'date', b'Date Field'), (b'href', b'Link to server or external'), (b'imghref', b'Image link to embed'), (b'decimal', b'Decimal field'), (b'boolean', b'checkbox'), (b'related', b'TEST'), (b'file', b'File Upload')]),
            preserve_default=True,
        ),
    ]
