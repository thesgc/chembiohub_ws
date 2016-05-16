# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.db import models, migrations
from rdkit import Chem

from rdkit.Chem import Descriptors as NewDescriptors

def recalculate_properties(apps, stuff):
    CompoundStructures = apps.get_model("chembl_business_model", "CompoundStructures")
    CompoundProperties = apps.get_model("chembl_business_model", "CompoundProperties")

    for struc in CompoundStructures.objects.all():
        molecule = struc.molecule
        if not hasattr(molecule, 'compoundproperties'):
            prop = CompoundProperties(molecule=molecule)
        else:
            prop = molecule.compoundproperties
        mol = Chem.MolFromMolBlock(str(struc.molfile))

        prop.full_mwt = NewDescriptors.MolWt(mol)
        prop.save()

class Migration(migrations.Migration):

    dependencies = [
        ('chembl_business_model', '0001_initial'),
    ]

    operations = [
        migrations.RunPython(recalculate_properties)
    ]
