# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0015_auto_20150914_0954'),
    ]

    operations = [
        migrations.AddField(
            model_name='datatype',
            name='uri',
            field=models.CharField(default=b'', max_length=1000),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='datatype',
            name='version',
            field=models.CharField(default=b'', max_length=10),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='field_type',
            field=models.CharField(default=b'char', max_length=15, choices=[(b'text', b'Short text field'), (b'char', b'Short text field'), (b'textarea', b'Full text'), (b'uiselect', b'Choice field'), (b'integer', b'Integer field'), (b'number', b'Decimal field'), (
                b'uiselecttag', b'Choice allowing create'), (b'uiselecttags', b'Tags field allowing create'), (b'percentage', b'Percentage field'), (b'date', b'Date Field'), (b'href', b'Link to server or external'), (b'imghref', b'Image link to embed')]),
            preserve_default=True,
        ),
    ]
