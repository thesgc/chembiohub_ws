# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations

#  alter table cbh_chembl_model_extension_project rename to cbh_core_model_project;
# ALTER TABLE
# cbh_chembl=# alter table cbh_chembl_model_extension_customfieldconfig rename to cbh_core_model_customfieldconfig;
# ALTER TABLE
# cbh_chembl=# alter table cbh_chembl_model_extension_pinnedcustomfield rename to cbh_core_model_pinnedcustomfield;
# ALTER TABLE
# cbh_chembl=# alter table cbh_chembl_model_extension_
# cbh_chembl_model_extension_cbhcompoundbatch          cbh_chembl_model_extension_cbhcompoundmultiplebatch  cbh_chembl_model_extension_projecttype               cbh_chembl_model_extension_skinningconfig
# cbh_chembl=# alter table cbh_chembl_model_extension_projecttype rename to cbh_core_model_projecttype;
# ALTER TABLE
# cbh_chembl=# alter table cbh_chembl_model_extension_skinningconfig
# rename to cbh_core_model_skinningconfig;


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0023_auto_20150806_0206'),
        ('cbh_core_model', '0001_initial'),

        ('chembl_core_model', '0007_auto_20150806_0115'),
    ]

    state_operations = [

        migrations.AlterField(
            model_name='cbhcompoundbatch',
            name='project',
            field=models.ForeignKey(
                default=None, blank=True, to='cbh_core_model.Project', null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='cbhcompoundmultiplebatch',
            name='project',
            field=models.ForeignKey(
                default=None, blank=True, to='cbh_core_model.Project', null=True),
            preserve_default=True,
        ),

    ]

    operations = [
        # By running only state operations, we are making Django think it has
        # applied this migration to the database. In reality, we renamed a
        # "cars_tires" table to "tires_tires" earlier.
        migrations.SeparateDatabaseAndState(state_operations=state_operations)
    ]
