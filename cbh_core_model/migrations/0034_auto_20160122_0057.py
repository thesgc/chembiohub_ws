# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


def add_skinning(apps, foo):
    """Adding the default singleton object for skinningconfig so it does not have to be added in the admin"""
    SkinningConfig = apps.get_model("cbh_core_model", "skinningconfig")
    SkinningConfig.objects.get_or_create(pk=1)


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0033_auto_20160121_1126'),
    ]

    operations = [
        migrations.RunPython(add_skinning)
    ]
