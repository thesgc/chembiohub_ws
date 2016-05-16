# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0033_auto_20160404_1128'),
    ]

    operations = [

        migrations.RunSQL(
            "alter sequence cbh_chembl_model_extension_project_id_seq rename to cbh_core_model_project_id_seq;"),
        migrations.RunSQL(
            "alter table cbh_chembl_model_extension_customfieldconfig_id_seq rename to cbh_core_model_customfieldconfig_id_seq;"),
        migrations.RunSQL(
            "alter table cbh_chembl_model_extension_pinnedcustomfield_id_seq rename to cbh_core_model_pinnedcustomfield_id_seq;"),
        migrations.RunSQL(
            "alter table cbh_chembl_model_extension_projecttype_id_seq rename to cbh_core_model_projecttype_id_seq;"),
        migrations.RunSQL(
            "alter table cbh_chembl_model_extension_skinningconfig_id_seq rename to cbh_core_model_skinningconfig_id_seq;"),
    ]
