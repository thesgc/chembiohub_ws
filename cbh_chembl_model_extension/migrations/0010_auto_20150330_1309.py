# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from rdkit.Chem import AllChem, MolFromMolBlock, MolToMolBlock


def calculate_coords(apps, schema_editor):
    # We can't import the Person model directly as it may be a newer
    # version than this migration expects. We use the historical version.
    Batch = apps.get_model("cbh_chembl_model_extension", "CBHCompoundBatch")
    for field in Batch.objects.all():
        mol = MolFromMolBlock(field.ctab)
        AllChem.Compute2DCoords(mol)
        try:
            field.ctab = MolToMolBlock(mol, includeStereo=True)
        except:
            print "test"
        field.save()


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0009_auto_20150330_1250'),
    ]

    operations = [
        migrations.RunPython(calculate_coords),
    ]
