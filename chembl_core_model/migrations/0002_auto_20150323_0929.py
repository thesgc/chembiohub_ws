# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import chembl_core_db.db.customFields


class Migration(migrations.Migration):

    dependencies = [
        ('chembl_core_model', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='CompoundStructuralAlerts',
            fields=[
                ('cpd_str_alert_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Primary key.', serialize=False, primary_key=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table':'compound_structural_alerts'
            },
        ),
        migrations.CreateModel(
            name='FracClassification',
            fields=[
                ('frac_class_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Unique numeric primary key for each level5 code', serialize=False, primary_key=True, blank=True)),
                ('active_ingredient', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of active ingredient (fungicide) classified by FRAC', max_length=500)),
                ('level1', chembl_core_db.db.customFields.ChemblCharField(help_text='Mechanism of action code assigned by FRAC', max_length=2)),
                ('level1_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of mechanism of action', max_length=2000)),
                ('level2', chembl_core_db.db.customFields.ChemblCharField(help_text='Target site code assigned by FRAC', max_length=2)),
                ('level2_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of target provided by FRAC', max_length=2000, null=True, blank=True)),
                ('level3', chembl_core_db.db.customFields.ChemblCharField(help_text='Group number assigned by FRAC', max_length=6)),
                ('level3_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of group provided by FRAC', max_length=2000, null=True, blank=True)),
                ('level4', chembl_core_db.db.customFields.ChemblCharField(help_text='Number denoting the chemical group (number not assigned by FRAC)', max_length=7)),
                ('level4_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Chemical group name provided by FRAC', max_length=2000, null=True, blank=True)),
                ('level5', chembl_core_db.db.customFields.ChemblCharField(help_text='A unique code assigned to each ingredient (based on the level 1-4 FRAC classification, but not assigned by IRAC)', unique=True, max_length=8)),
                ('frac_code', chembl_core_db.db.customFields.ChemblCharField(help_text='The official FRAC classification code for the ingredient', max_length=4)),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'frac_classification'
            },
        ),
        migrations.CreateModel(
            name='HracClassification',
            fields=[
                ('hrac_class_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Unique numeric primary key for each level3 code', serialize=False, primary_key=True, blank=True)),
                ('active_ingredient', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of active ingredient (herbicide) classified by HRAC', max_length=500)),
                ('level1', chembl_core_db.db.customFields.ChemblCharField(help_text='HRAC group code - denoting mechanism of action of herbicide', max_length=2)),
                ('level1_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of mechanism of action provided by HRAC', max_length=2000)),
                ('level2', chembl_core_db.db.customFields.ChemblCharField(help_text='Indicates a chemical family within a particular HRAC group (number not assigned by HRAC)', max_length=3)),
                ('level2_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of chemical family provided by HRAC', max_length=2000, null=True, blank=True)),
                ('level3', chembl_core_db.db.customFields.ChemblCharField(help_text='A unique code assigned to each ingredient (based on the level 1 and 2 HRAC classification, but not assigned by HRAC)', unique=True, max_length=5)),
                ('hrac_code', chembl_core_db.db.customFields.ChemblCharField(help_text='The official HRAC classification code for the ingredient', max_length=2)),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'hrac_classification'
            },
        ),
        migrations.CreateModel(
            name='IracClassification',
            fields=[
                ('irac_class_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Unique numeric primary key for each level4 code', serialize=False, primary_key=True, blank=True)),
                ('active_ingredient', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of active ingredient (insecticide) classified by IRAC', max_length=500)),
                ('level1', chembl_core_db.db.customFields.ChemblCharField(help_text='Class of action e.g., nerve action, energy metabolism (code not assigned by IRAC)', max_length=1, choices=[(b'A', b'A'), (b'B', b'B'), (b'C', b'C'), (b'D', b'D'), (b'E', b'E'), (b'M', b'M'), (b'U', b'U')])),
                ('level1_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of class of action, as provided by IRAC', max_length=2000, choices=[(b'ENERGY METABOLISM', b'ENERGY METABOLISM'), (b'GROWTH REGULATION', b'GROWTH REGULATION'), (b'LIPID SYNTHESIS, GROWTH REGULATION', b'LIPID SYNTHESIS, GROWTH REGULATION'), (b'MISCELLANEOUS', b'MISCELLANEOUS'), (b'NERVE ACTION', b'NERVE ACTION'), (b'NERVE AND MUSCLE ACTION', b'NERVE AND MUSCLE ACTION'), (b'UNKNOWN', b'UNKNOWN')])),
                ('level2', chembl_core_db.db.customFields.ChemblCharField(help_text='IRAC main group code denoting primary site/mechanism of action', max_length=3)),
                ('level2_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of site/mechanism of action provided by IRAC', max_length=2000)),
                ('level3', chembl_core_db.db.customFields.ChemblCharField(help_text='IRAC sub-group code denoting chemical class of insecticide', max_length=6)),
                ('level3_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of chemical class or exemplifying ingredient provided by IRAC', max_length=2000)),
                ('level4', chembl_core_db.db.customFields.ChemblCharField(help_text='A unique code assigned to each ingredient (based on the level 1, 2 and 3 IRAC classification, but not assigned by IRAC)', unique=True, max_length=8)),
                ('irac_code', chembl_core_db.db.customFields.ChemblCharField(help_text='The official IRAC classification code for the ingredient', max_length=3)),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table':'irac_classification'
            },
        ),
        migrations.CreateModel(
            name='MoleculeFracClassification',
            fields=[
                ('mol_frac_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
                ('frac_class', models.ForeignKey(help_text='Foreign key to frac_classification table showing the mechanism of action classification of the compound.', to='chembl_core_model.FracClassification')),
                ('molecule', models.ForeignKey(db_column=b'molregno', to='chembl_core_model.MoleculeDictionary', help_text='Foreign key to molecule_dictionary, showing the compound to which the classification applies.')),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'molecule_frac_classification'
            },
        ),
        migrations.CreateModel(
            name='MoleculeHracClassification',
            fields=[
                ('mol_hrac_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key', serialize=False, primary_key=True, blank=True)),
                ('hrac_class', models.ForeignKey(help_text='Foreign key to hrac_classification table showing the classification for the compound.', to='chembl_core_model.HracClassification')),
                ('molecule', models.ForeignKey(db_column=b'molregno', to='chembl_core_model.MoleculeDictionary', help_text='Foreign key to molecule_dictionary, showing the compound to which this classification applies.')),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'molecule_hrac_classification'
            },
        ),
        migrations.CreateModel(
            name='MoleculeIracClassification',
            fields=[
                ('mol_irac_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
                ('irac_class', models.ForeignKey(help_text='Foreign key to the irac_classification table showing the mechanism of action classification for the compound.', to='chembl_core_model.IracClassification')),
                ('molecule', models.ForeignKey(db_column=b'molregno', to='chembl_core_model.MoleculeDictionary', help_text='Foreign key to the molecule_dictionary table, showing the compound to which the classification applies.')),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'molecule_irac_classification'
            },
        ),
        migrations.CreateModel(
            name='PatentUseCodes',
            fields=[
                ('patent_use_code', chembl_core_db.db.customFields.ChemblCharField(help_text='Primary key. Patent use code from FDA Orange Book', max_length=8, serialize=False, primary_key=True)),
                ('definition', chembl_core_db.db.customFields.ChemblCharField(help_text='Definition for the patent use code, from FDA Orange Book.', max_length=500)),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'patent_use_codes'
            },
        ),
        migrations.CreateModel(
            name='ProductPatents',
            fields=[
                ('prod_pat_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key', serialize=False, primary_key=True, blank=True)),
                ('patent_no', chembl_core_db.db.customFields.ChemblCharField(help_text='Patent numbers as submitted by the applicant holder for patents covered by the statutory provisions', max_length=11)),
                ('patent_expire_date', chembl_core_db.db.customFields.ChemblDateField(help_text='Date the patent expires as submitted by the applicant holder including applicable extensions')),
                ('drug_substance_flag', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Patents submitted on FDA Form 3542 and listed after August 18, 2003 may have a drug substance flag set to 1, indicating the sponsor submitted the patent as claiming the drug substance')),
                ('drug_product_flag', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Patents submitted on FDA Form 3542 and listed after August 18, 2003 may have a drug product flag set to 1, indicating the sponsor submitted the patent as claiming the drug product')),
                ('delist_flag', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Sponsor has requested patent be delisted if set to 1.  This patent has remained listed because, under Section 505(j)(5)(D)(i) of the Act, a first applicant may retain eligibility for 180-day exclusivity based on a paragraph IV certification to this patent for a certain period.  Applicants under Section 505(b)(2) are not required to certify to patents where this flag is set to 1')),
                ('in_products', chembl_core_db.db.customFields.ChemblPositiveIntegerField(default=0, help_text='Indicates whether the PRODUCT_ID can be found in the PRODUCTS table (where set to 1)')),
                ('patent_use_code', models.ForeignKey(db_column=b'patent_use_code', blank=True, to='chembl_core_model.PatentUseCodes', help_text='Code to designate a use patent that covers the approved indication or use of a drug product', null=True)),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table':'product_patents'
            },
        ),
        migrations.CreateModel(
            name='StructuralAlerts',
            fields=[
                ('alert_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Primary key. Unique identifier for the structural alert', serialize=False, primary_key=True)),
                ('alert_name', chembl_core_db.db.customFields.ChemblCharField(help_text='A name for the structural alert', max_length=100)),
                ('smarts', chembl_core_db.db.customFields.ChemblCharField(help_text='SMARTS defining the structural feature that is considered to be an alert', max_length=4000)),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'structural_alerts'
            },
        ),
        migrations.CreateModel(
            name='StructuralAlertSets',
            fields=[
                ('alert_set_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Unique ID for the structural alert set', serialize=False, primary_key=True, choices=[(1, b'1'), (2, b'2'), (3, b'3'), (4, b'4'), (5, b'5'), (6, b'6'), (7, b'7'), (8, b'8')])),
                ('set_name', chembl_core_db.db.customFields.ChemblCharField(help_text='Name (or origin) of the structural alert set', unique=True, max_length=100, choices=[(b'BMS', b'BMS'), (b'Dundee', b'Dundee'), (b'Glaxo', b'Glaxo'), (b'Inpharmatica', b'Inpharmatica'), (b'LINT', b'LINT'), (b'MLSMR', b'MLSMR'), (b'PAINS', b'PAINS'), (b'SureChEMBL', b'SureChEMBL')])),
                ('priority', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Priority assigned to the structural alert set for display on the ChEMBL interface (priorities >=4 are shown by default).', choices=[(1, b'1'), (2, b'2'), (3, b'3'), (4, b'4'), (5, b'5'), (6, b'6'), (7, b'7'), (8, b'8')])),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'structural_alert_sets'
            },
        ),
        migrations.RemoveField(
            model_name='activities',
            name='activity_type',
        ),
        migrations.RemoveField(
            model_name='compoundrecords',
            name='old_compound_key',
        ),
        migrations.RemoveField(
            model_name='drugmechanism',
            name='uniprot_accessions',
        ),
        migrations.RemoveField(
            model_name='targetdictionary',
            name='in_starlite',
        ),
        migrations.RemoveField(
            model_name='targetdictionary',
            name='popularity',
        ),
        migrations.AddField(
            model_name='biotherapeutics',
            name='helm_notation',
            field=chembl_core_db.db.customFields.ChemblCharField(help_text='Sequence notation generated according to the HELM standard (http://www.openhelm.org/home). Currently for peptides only', max_length=4000, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='celldictionary',
            name='chembl',
            field=models.ForeignKey(blank=True, to='chembl_core_model.ChemblIdLookup', help_text='ChEMBL identifier for the cell (used in web interface etc)', null=True),
        ),
        migrations.AddField(
            model_name='celldictionary',
            name='cl_lincs_id',
            field=chembl_core_db.db.customFields.ChemblCharField(help_text='Cell ID used in LINCS (Library of Integrated Network-based Cellular Signatures)', max_length=8, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='compoundproperties',
            name='hba_lipinski',
            field=chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text="Number of hydrogen bond acceptors calculated according to Lipinski's original rules (i.e., N + O count))", null=True, blank=True),
        ),
        migrations.AddField(
            model_name='compoundproperties',
            name='hbd_lipinski',
            field=chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text="Number of hydrogen bond donors calculated according to Lipinski's original rules (i.e., NH + OH count)", null=True, blank=True),
        ),
        migrations.AddField(
            model_name='compoundproperties',
            name='num_lipinski_ro5_violations',
            field=chembl_core_db.db.customFields.ChemblPositiveIntegerField(blank=True, help_text="Number of violations of Lipinski's rule of five using HBA_LIPINSKI and HBD_LIPINSKI counts", null=True, choices=[(0, b'0'), (1, b'1'), (2, b'2'), (3, b'3'), (4, b'4')]),
        ),
        #migrations.AlterField(
        #    model_name='celldictionary',
        #    name='downgraded',
        #    field=chembl_core_db.db.customFields.ChemblNullableBooleanField(default=False, help_text='Indicates the cell line has been removed (if set to 1)'),
        #),
        migrations.AlterField(
            model_name='chemblidlookup',
            name='entity_id',
            field=chembl_core_db.db.customFields.ChemblIntegerField(help_text='Primary key for that entity in corresponding table (e.g., molregno for compounds, tid for targets)'),
        ),
        migrations.AlterField(
            model_name='chemblidlookup',
            name='entity_type',
            field=chembl_core_db.db.customFields.ChemblCharField(help_text='Type of entity (e.g., COMPOUND, ASSAY, TARGET)', max_length=50,)
        ),
        migrations.AlterField(
            model_name='chemblidlookup',
            name='status',
            field=chembl_core_db.db.customFields.ChemblCharField(default='ACTIVE', help_text='Indicates whether the status of the entity within the database - ACTIVE, INACTIVE (downgraded), OBS (obsolete/removed).', max_length=10, choices=[(b'ACTIVE', b'ACTIVE'), (b'INACTIVE', b'INACTIVE'), (b'OBS', b'OBS')]),
        ),
        #migrations.AlterField(
        #    model_name='compoundproperties',
        #    name='acd_most_bpka',
        #    field=chembl_core_db.db.customFields.ChemblPositiveDecimalField(help_text='The most basic pKa calculated using ACDlabs v12.01', null=True, max_digits=9, decimal_places=2, blank=True),
       # ),
        migrations.AlterField(
            model_name='compoundproperties',
            name='med_chem_friendly',
            field=chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=3, null=True, help_text='DEPRECATED. Replaced by new structural alerts tables. Will be removed in future releases.', choices=[(b'Y', b'Yes'), (b'N', b'No')]),
        ),
        migrations.AlterField(
            model_name='compoundproperties',
            name='num_alerts',
            field=chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Number of structural alerts used for QED calculation (as defined by Brenk et al., ChemMedChem 2008)', null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='compoundproperties',
            name='num_ro5_violations',
            field=chembl_core_db.db.customFields.ChemblPositiveIntegerField(blank=True, help_text="Number of violations of Lipinski's rule-of-five, using HBA and HBD definitions", null=True, db_index=True, choices=[(0, b'0'), (1, b'1'), (2, b'2'), (3, b'3'), (4, b'4')]),
        ),
        #migrations.AlterField(
        #    model_name='defineddailydose',
        #    name='ddd_value',
        #    field=chembl_core_db.db.customFields.ChemblPositiveDecimalField(help_text='Value of defined daily dose', null=True, max_digits=9, decimal_places=2, blank=True),
        #),
        migrations.AlterField(
            model_name='drugmechanism',
            name='selectivity_comment',
            field=chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=100, null=True, help_text='Additional comments regarding the selectivity of the drug', ),
        ),
        migrations.AlterField(
            model_name='products',
            name='information_source',
            field=chembl_core_db.db.customFields.ChemblCharField(help_text='Source of the product information (e.g., Orange Book)', max_length=100, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='products',
            name='product_class',
            field=chembl_core_db.db.customFields.ChemblCharField(max_length=30, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='proteinclassification',
            name='class_level',
            field=chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Level of the class within the hierarchy (level 1 = top level classification)', choices=[(0, b'0'), (1, b'1'), (2, b'2'), (3, b'3'), (4, b'4'), (5, b'5'), (6, b'6'), (7, b'7'), (8, b'8')]),
        ),
       # migrations.AlterField(
       #     model_name='proteinclassification',
       #     name='downgraded',
       #     field=chembl_core_db.db.customFields.ChemblBooleanField(default=False),
       # ),
        #migrations.AlterField(
        #    model_name='targetcomponents',
        #    name='homologue',
        #    field=chembl_core_db.db.customFields.ChemblPositiveIntegerField(default=0, help_text='Indicates that the given component is a homologue of the correct component (e.g., from a different species) when set to 1. This may be the case if the sequence for the correct protein/nucleic acid cannot be found in sequence databases. A value of 2 indicates that the sequence given is a representative of a species group, e.g., an E. coli protein to represent the target of a broad-spectrum antibiotic.', ),
        #),
        migrations.AlterField(
            model_name='targetcomponents',
            name='relationship',
            field=chembl_core_db.db.customFields.ChemblCharField( max_length=20, ),
        ),
        #migrations.AlterField(
        #    model_name='targetdictionary',
        #    name='downgraded',
        #    field=chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Flag to indicate that the target is downgraded (if equal to 1)'),
        #),
        migrations.AlterField(
            model_name='targetdictionary',
            name='pref_name',
            field=chembl_core_db.db.customFields.ChemblCharField(help_text='Preferred target name: manually curated', max_length=200, db_index=True),
        ),
       # migrations.AlterField(
       #     model_name='targetdictionary',
       #     name='species_group_flag',
       #     field=chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text="Flag to indicate whether the target represents a group of species, rather than an individual species (e.g., 'Bacterial DHFR'). Where set to 1, indicates that any associated target components will be a representative, rather than a comprehensive set."),
       # ),
        migrations.AlterField(
            model_name='targetrelations',
            name='related_target',
            field=models.ForeignKey(related_name='from_target', db_column=b'related_tid', to='chembl_core_model.TargetDictionary', help_text='Identifier for the target that is related to the target of interest (foreign key to target_dicitionary table)'),
        ),
        migrations.AddField(
            model_name='structuralalerts',
            name='alert_set',
            field=models.ForeignKey(choices=[(1, b'1'), (2, b'2'), (3, b'3'), (4, b'4'), (5, b'5'), (6, b'6'), (7, b'7'), (8, b'8')], to='chembl_core_model.StructuralAlertSets', help_text='Foreign key to structural_alert_sets table indicating which set this particular alert comes from'),
        ),
        migrations.AddField(
            model_name='productpatents',
            name='product',
            field=models.ForeignKey(help_text='Foreign key to products table - FDA application number for the product', to='chembl_core_model.Products'),
        ),
        migrations.AddField(
            model_name='compoundstructuralalerts',
            name='alert',
            field=models.ForeignKey(help_text='Foreign key to the structural_alerts table. The particular alert that has been identified in this compound.', to='chembl_core_model.StructuralAlerts'),
        ),
        migrations.AddField(
            model_name='compoundstructuralalerts',
            name='molecule',
            field=models.ForeignKey(db_column=b'molregno', to='chembl_core_model.MoleculeDictionary', help_text='Foreign key to the molecule_dictionary. The compound for which the structural alert has been found.'),
        ),
        migrations.AlterUniqueTogether(
            name='productpatents',
            unique_together=set([('product', 'patent_no', 'patent_expire_date', 'patent_use_code')]),
        ),
        migrations.AlterUniqueTogether(
            name='moleculeiracclassification',
            unique_together=set([('irac_class', 'molecule')]),
        ),
        migrations.AlterUniqueTogether(
            name='moleculehracclassification',
            unique_together=set([('hrac_class', 'molecule')]),
        ),
        migrations.AlterUniqueTogether(
            name='moleculefracclassification',
            unique_together=set([('frac_class', 'molecule')]),
        ),
        migrations.AlterUniqueTogether(
            name='compoundstructuralalerts',
            unique_together=set([('molecule', 'alert')]),
        ),
    ]
