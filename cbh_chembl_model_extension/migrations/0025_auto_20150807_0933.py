# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0024_auto_20150806_0115'),
    ]

    state_operations = [
        migrations.RemoveField(
            model_name='customfieldconfig',
            name='created_by',
        ),
        migrations.RemoveField(
            model_name='pinnedcustomfield',
            name='custom_field_config',
        ),
        migrations.DeleteModel(
            name='PinnedCustomField',
        ),
        migrations.RemoveField(
            model_name='project',
            name='created_by',
        ),
        migrations.RemoveField(
            model_name='project',
            name='custom_field_config',
        ),
        migrations.DeleteModel(
            name='CustomFieldConfig',
        ),
        migrations.RemoveField(
            model_name='project',
            name='project_type',
        ),
        migrations.DeleteModel(
            name='Project',
        ),
        migrations.DeleteModel(
            name='ProjectType',
        ),
        migrations.DeleteModel(
            name='SkinningConfig',
        ),
    ]

    operations = [
        # By running only state operations, we are making Django think it has
        # applied this migration to the database. In reality, we renamed a
        # "cars_tires" table to "tires_tires" earlier.
        migrations.SeparateDatabaseAndState(state_operations=state_operations)
    ]
