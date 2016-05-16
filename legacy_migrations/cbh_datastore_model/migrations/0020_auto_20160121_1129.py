# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model', '0019_auto_20160121_1126'),
    ]

    operations = [
        migrations.AlterField(
            model_name='attachment',
            name='attachment_custom_field_config',
            field=models.ForeignKey(default=None, blank=True, to='cbh_core_model.CustomFieldConfig', help_text=b'The schema of the table in this attachment', null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='attachment',
            name='chosen_data_form_config',
            field=models.ForeignKey(default=None, blank=True, to='cbh_core_model.DataFormConfig', help_text=b'The template data form config whose last level corresponds to the data being added', null=True),
            preserve_default=True,
        ),
    ]
