# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.template.defaultfilters import slugify


def create_keys(apps, schema_editor):
    # We can't import the Person model directly as it may be a newer
    # version than this migration expects. We use the historical version.
    PinnedCustomField = apps.get_model(
        "cbh_chembl_model_extension", "PinnedCustomField")
    for field in PinnedCustomField.objects.all():
        field.field_key = slugify(field.name).replace("-", "_")
        field.save()


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0005_auto_20150330_1047'),
    ]

    operations = [
        migrations.RunPython(create_keys),
    ]
