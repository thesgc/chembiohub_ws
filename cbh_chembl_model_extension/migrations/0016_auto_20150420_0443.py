# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


def migrate_multiple_batch_data(apps, schema_editor):
    # We can't import the Person model directly as it may be a newer
    # version than this migration expects. We use the historical version.
    Batch = apps.get_model("cbh_chembl_model_extension", "CBHCompoundBatch")
    MultipleBatch = apps.get_model(
        "cbh_chembl_model_extension", "CBHCompoundMultipleBatch")

    for batch in Batch.objects.all():
        if batch.multiple_batch_id:

            mb = MultipleBatch.objects.get(pk=batch.multiple_batch_id)
            mb.created_by = batch.created_by
            mb.project = batch.project
            mb.save()


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0015_auto_20150417_1631'),
    ]

    operations = [
        migrations.RunPython(migrate_multiple_batch_data),

    ]
