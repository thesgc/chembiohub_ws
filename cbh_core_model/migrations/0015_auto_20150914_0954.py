# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0014_dataformconfig_parent'),
    ]

    operations = [
        migrations.AddField(
            model_name='pinnedcustomfield',
            name='pinned_for_datatype',
            field=models.ForeignKey(
                default=None, blank=True, to='cbh_core_model.DataType', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='pinnedcustomfield',
            name='standardised_alias',
            field=models.ForeignKey(related_name='alias_mapped_from', default=None,
                                    blank=True, to='cbh_core_model.PinnedCustomField', null=True),
            preserve_default=True,
        ),
        
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='custom_field_config',
            field=models.ForeignKey(related_name='pinned_custom_field', default=None,
                                    blank=True, to='cbh_core_model.CustomFieldConfig', null=True),
            preserve_default=True,
        ),
    ]
