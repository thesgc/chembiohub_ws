# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import chembl_core_db.db.customFields


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0023_auto_20150806_0206'),

        ('cbh_core_model', '0001_initial'),

        ('chembl_core_model', '0006_auto_20150521_0906'),
    ]

    state_operations = [

        migrations.AlterField(
            model_name='moleculedictionary',
            name='project',
            field=models.ForeignKey(blank=True, to='cbh_core_model.Project', null=True),
            preserve_default=True,
        ),
     
    ]


    operations = [

        # By running only state operations, we are making Django think it has
        # applied this migration to the database. In reality, we renamed a
        # "cars_tires" table to "tires_tires" earlier.
        migrations.SeparateDatabaseAndState(state_operations=state_operations)
    ]