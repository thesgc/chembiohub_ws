# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import chembl_core_db.db.customFields


class Migration(migrations.Migration):

    dependencies = [
        ('chembl_core_model', '0005_auto_20150323_1105'),
    ]

    operations = [
        # migrations.AlterField(
        #     model_name='chemblidlookup',
        #     name='entity_id',
        #     field=chembl_core_db.db.customFields.ChemblIntegerField(default=-1, help_text='Primary key for that entity in corresponding table (e.g., molregno for compounds, tid for targets)', null=True),
        #     preserve_default=True,
        # ),
        # migrations.AlterField(
        #     model_name='chemblidlookup',
        #     name='entity_type',
        #     field=chembl_core_db.db.customFields.ChemblCharField(help_text='Type of entity (e.g., COMPOUND, ASSAY, TARGET)', max_length=50, choices=[(b'ASSAY', b'ASSAY'), (b'CELL', b'CELL'), (b'COMPOUND', b'COMPOUND'), (b'DOCUMENT', b'DOCUMENT'), (b'TARGET', b'TARGET'), (b'BLINDED_BATCH', b'BLINDED_BATCH')]),
        #     preserve_default=True,
        # ),
        # migrations.AlterField(
        #     model_name='moleculedictionary',
        #     name='forced_reg_index',
        #     field=chembl_core_db.db.customFields.ChemblPositiveIntegerField(default=0, help_text='Number of times this structure key has been forced to be registered', db_index=True),
        #     preserve_default=True,
        # ),
        # migrations.AlterField(
        #     model_name='moleculedictionary',
        #     name='public',
        #     field=chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Whether this molecule has been marked as public'),
        #     preserve_default=True,
        # ),
        # migrations.AlterField(
        #     model_name='moleculedictionary',
        #     name='structure_type',
        #     field=chembl_core_db.db.customFields.ChemblCharField(default=b'MOL', help_text='Indications whether the molecule has a small molecule structure or a protein sequence (MOL indicates an entry in the compound_structures table, SEQ indications an entry in the protein_therapeutics table, NONE indicates an entry in neither table, e.g., structure unknown)', max_length=10, choices=[(b'NONE', b'NONE'), (b'MOL', b'MOL'), (b'SEQ', b'SEQ'), (b'BOTH', b'BOTH'), (b'BLIND_ID', b'Blind ID')]),
        #     preserve_default=True,
        # ),
    ]
