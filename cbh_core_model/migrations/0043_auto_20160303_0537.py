# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations

def update_field_type(apps, blah):
    pcf = apps.get_model('cbh_core_model.PinnedCustomField')
    for p in pcf.objects.filter(field_type="object"):
        p.field_type = "file"
        p.save()

class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0042_auto_20160303_0528'),
    ]

    operations = [
    	migrations.RunPython(update_field_type)
    ]
