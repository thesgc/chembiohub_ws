# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations

def migrate_perms(apps, schema_editor):
    ContentType = apps.get_model("contenttypes", "contenttype")
    Permission = apps.get_model("auth", "Permission")
    ct, created = ContentType.objects.get_or_create(
            app_label="_can_see", model=ContentType)
      
    pm = Permission.objects.get_or_create( content_type_id=ct.id, codename="no_chemreg")
    pm = Permission.objects.get_or_create( content_type_id=ct.id, codename="no_assayreg")


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0022_auto_20150925_0337'),
    ]

    operations = [
        migrations.RunPython(migrate_perms)
    ]
