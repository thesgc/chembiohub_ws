# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations

from django.core.exceptions import ObjectDoesNotExist




def migrate_old_permissions_to_new_ones(apps, schema_editor):
    """Idempotent function to migrate the old permissions to the new ones. Old custom contentypes are not deleted for now
    This will allow the app to be migrated to django 1.9 and keeps the permissions system in line with django's meaning 
    functions such as user.is_administrator dontate all permissions to the user properly

    """

    Permission = apps.get_model("auth", "Permission")
    Project = apps.get_model("cbh_core_model", "Project")
    ContentType = apps.get_model("contenttypes", "ContentType")
    try:
        new_ct = ContentType.objects.get(app_label="cbh_core_model", model="project")
        skin_ct = ContentType.objects.get(app_label="cbh_core_model", model="skinningconfig")
        from cbh_core_model.models import get_permission_name, get_permission_codename
    
        for perm in Permission.objects.all():
            ct = perm.content_type
            if ct.app_label.isdigit():
                try:
                    project = Project.objects.get(id=int(ct.app_label))
                    if perm.codename == "admin":
                        perm.codename = "owner"
                    perm.name = get_permission_name(project.name, perm.codename)
                    perm.codename = get_permission_codename(project.id, perm.codename)
                    
                    perm.content_type = new_ct
                    perm.content_type_id = new_ct.id
                    perm.save()
                except ObjectDoesNotExist:
                    pass
            if ct.app_label == "_can_see":
                perm.content_type = skin_ct
                perm.content_type_id = skin_ct.id
                perm.save()#
    except ObjectDoesNotExist:
        print("No contenttypes therefore nothing to migrate")





class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0030_auto_20151215_0548'),
    ]

    operations = [
        migrations.RunPython(migrate_old_permissions_to_new_ones,)
    ]
