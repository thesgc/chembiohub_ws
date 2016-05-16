# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    def create_default_datapoint(apps, b):
        User = apps.get_model("auth.User")
        User.objects.get_or_create(username=-1, id=-1)


        CFC = apps.get_model("cbh_core_model.CustomFieldConfig")
        CFC.objects.get_or_create(
                id=-1, created_by_id=-1, name="-1 default do not delete")


        DataPoint = apps.get_model("cbh_datastore_model.DataPoint")
        dp = DataPoint.objects.get_or_create(
            custom_field_config_id=-1, created_by_id=-1, project_data={}, supplementary_data={})


        DataPointClassification = apps.get_model(
                "cbh_datastore_model.DataPointClassification")
        DataPointClassification.objects.all().delete()


    dependencies = [
        ('cbh_datastore_model', '0010_auto_20150810_0917'),
    ]

    operations = [
        migrations.RunPython(create_default_datapoint),
    ]
