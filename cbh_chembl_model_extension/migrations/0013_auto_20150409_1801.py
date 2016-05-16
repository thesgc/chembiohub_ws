# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from copy import deepcopy


def migrate_unknown_custom_fields(apps, schema_editor):
    # We can't import the Person model directly as it may be a newer
    # version than this migration expects. We use the historical version.
    Batch = apps.get_model("cbh_chembl_model_extension", "CBHCompoundBatch")

    for batch in Batch.objects.all():
        batch.uncurated_fields = deepcopy(batch.custom_fields)
        batch.custom_fields = {}
        batch.save()


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0012_auto_20150408_0740'),
    ]

    operations = [
        migrations.RenameField(
            model_name='cbhcompoundbatch',
            old_name='viewable_by',
            new_name='uncurated_fields',
        ),
        migrations.RunPython(migrate_unknown_custom_fields),

    ]
