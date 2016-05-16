# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0022_auto_20150721_1031'),
    ]

    operations = [
        migrations.RunSQL(
            "alter table cbh_chembl_model_extension_project rename to cbh_core_model_project;"),
        migrations.RunSQL(
            "alter table cbh_chembl_model_extension_customfieldconfig rename to cbh_core_model_customfieldconfig;"),
        migrations.RunSQL(
            "alter table cbh_chembl_model_extension_pinnedcustomfield rename to cbh_core_model_pinnedcustomfield;"),
        migrations.RunSQL(
            "alter table cbh_chembl_model_extension_projecttype rename to cbh_core_model_projecttype;"),
        migrations.RunSQL(
            "alter table cbh_chembl_model_extension_skinningconfig rename to cbh_core_model_skinningconfig;"),


    ]
