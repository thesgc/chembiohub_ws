# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings
import chembl_core_db.db.customFields


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('cbh_chembl_model_extension', '__first__'),
        ('chembl_core_model', '0003_auto_20150323_1322'),
    ]

    operations = [
        migrations.AddField(
            model_name='moleculedictionary',
            name='created_by',
            field=models.ForeignKey(blank=True, to=settings.AUTH_USER_MODEL, null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='moleculedictionary',
            name='forced_reg_index',
            #field=chembl_core_db.db.customFields.ChemblPositiveIntegerField(default=0, help_text='Number of times this structure key has been forced to be registered', db_index=True),
            field=models.PositiveIntegerField(default=0, help_text='Number of times this structure key has been forced to be registered', db_index=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='moleculedictionary',
            name='forced_reg_reason',
            field=chembl_core_db.db.customFields.ChemblCharField(help_text='Reason for forced registration (e.g., known to be a stereoisomer)', max_length=200, null=True, blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='moleculedictionary',
            name='project',
            field=models.ForeignKey(blank=True, to='cbh_chembl_model_extension.Project', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='moleculedictionary',
            name='public',
            field=chembl_core_db.db.customFields.ChemblBooleanField(help_text='Whether this molecule has been marked as public'),
            #field=models.BooleanField(default=False, help_text='Whether this molecule has been marked as public'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='compoundstructures',
            name='standard_inchi',
            field=chembl_core_db.db.customFields.ChemblCharField(help_text='IUPAC standard InChI for the compound', max_length=4000, null=True, db_index=True, blank=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='moleculedictionary',
            name='molregno',
            field=models.AutoField(help_text='Internal Primary Key for the molecule', serialize=False, primary_key=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='moleculedictionary',
            name='structure_key',
            field=chembl_core_db.db.customFields.ChemblCharField(help_text='Note unique key removed for Chembiohub - Unique key for the structure/sequence (e.g., inchi_key or sequence md5sum) to help enforce non-redundancy.', max_length=27, null=True, db_index=True, blank=True),
            preserve_default=True,
        ),
        # migrations.AlterField(
        #     model_name='moleculedictionary',
        #     name='structure_type',
        #     field=chembl_core_db.db.customFields.ChemblCharField(default=b'MOL', help_text='Indications whether the molecule has a small molecule structure or a protein sequence (MOL indicates an entry in the compound_structures table, SEQ indications an entry in the protein_therapeutics table, NONE indicates an entry in neither table, e.g., structure unknown)', max_length=10, choices=[(b'NONE', b'NONE'), (b'MOL', b'MOL'), (b'SEQ', b'SEQ'), (b'BOTH', b'BOTH'), (b'BLIND_ID', b'Blind ID')]),
        #     preserve_default=True,
        # ),
        migrations.AlterUniqueTogether(
            name='moleculedictionary',
            unique_together=set([('structure_key', 'project', 'structure_type', 'forced_reg_index')]),
        ),
    ]
