# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


# def clone_field_list(field_list, custom_field_config_id):
#     for f in field_list:
#         f.id = None
#         f.pk = None
#         f.custom_field_config_id = custom_field_config_id
#         f.save()
#     return field_list

# def clone_custom_field_configs_by_project(apps, schema_editor):
#     '''If there is more than 1 project per custom field config then clone all the other custom field configs'''
#     CustomFieldConfig = apps.get_model("cbh_core_model", "CustomFieldConfig")
#     old_ids = []
#     for c in CustomFieldConfig.objects.all():
#         old_id = c.id
#         count_projects = c.project.count()
#         if count_projects > 1:
#             project_list = list(c.project.all())
#             field_list = list(c.pinned_custom_field.all())

#             for p in project_list:

#                 c.pk = None
#                 c.name = "%d__%s__project__config" % (p.id, p.name)
#                 c.project = []
#                 c.save()
#                 clone_field_list(field_list, c.id)

                
#                 p.custom_field_config_id = c.id

#                 p.save()
#             old_ids.append(old_id)





class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0031_auto_20160112_0307'),
    ]

    operations = [
        # migrations.RunPython(clone_custom_field_configs_by_project)
    ]
