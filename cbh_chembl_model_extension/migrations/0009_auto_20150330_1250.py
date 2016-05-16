# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0008_auto_20150330_1127'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='project',
            name='custom_field_config',
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='field_type',
            field=models.CharField(default=b'char', max_length=15, choices=[(b'checkboxes', b'Checkbox Fields'), (b'uiselect', b'Choice field'), (b'textarea', b'Full text'), (
                b'integer', b'Integer field'), (b'text', b'Short text field'), (b'uiselecttags', b'Tags field allowing create'), (b'number', b'Decimal field'), (b'uiselecttag', b'Choice allowing create')]),
            preserve_default=True,
        ),
    ]
