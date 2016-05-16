# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
import chembl_core_db.db.customFields





class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='ActionType',
            fields=[
                ('action_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Primary key. Type of action of the drug e.g., agonist, antagonist', max_length=50, serialize=False, primary_key=True)),
                ('description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of how the action type is used', max_length=200)),
                ('parent_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Higher-level grouping of action types e.g., positive vs negative action', max_length=50, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table' : 'action_type',                
            },
            bases=( models.Model,),
        ),
        migrations.CreateModel(
            name='Activities',
            fields=[
                ('activity_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Unique ID for the activity row', serialize=False, primary_key=True, blank=True)),
                ('activity_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Deprecated. Use published_activity_type or standard_type columns', max_length=250, null=True, db_index=True, blank=True)),
                ('standard_relation', chembl_core_db.db.customFields.ChemblCharField(choices=[(b'>', b'>'), (b'<', b'<'), (b'=', b'='), (b'~', b'~'), (b'<=', b'<='), (b'>=', b'>='), (b'<<', b'<<'), (b'>>', b'>>')], max_length=50, blank=True, help_text='Symbol constraining the activity value (e.g. >, <, =)', null=True, db_index=True)),
                ('published_value', chembl_core_db.db.customFields.ChemblNoLimitDecimalField(decimal_places=19, max_digits=38, blank=True, help_text='Datapoint value as it appears in the original publication.', null=True, db_index=True)),
                ('published_units', chembl_core_db.db.customFields.ChemblCharField(help_text='Units of measurement as they appear in the original publication', max_length=100, null=True, db_index=True, blank=True)),
                ('standard_value', chembl_core_db.db.customFields.ChemblNoLimitDecimalField(decimal_places=19, max_digits=38, blank=True, help_text='Same as PUBLISHED_VALUE but transformed to common units: e.g. mM concentrations converted to nM.', null=True, db_index=True)),
                ('standard_units', chembl_core_db.db.customFields.ChemblCharField(help_text="Selected 'Standard' units for data type: e.g. concentrations are in nM.", max_length=100, null=True, db_index=True, blank=True)),
                ('standard_flag', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Shows whether the standardised columns have been curated/set (1) or just default to the published data (0).')),
                ('standard_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Standardised version of the published_activity_type (e.g. IC50 rather than Ic-50/Ic50/ic50/ic-50)', max_length=250, null=True, db_index=True, blank=True)),
                ('updated_by', chembl_core_db.db.customFields.ChemblCharField(max_length=100, null=True, blank=True)),
                ('updated_on', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
                ('activity_comment', chembl_core_db.db.customFields.ChemblCharField(help_text="Describes non-numeric activities i.e. 'Slighty active', 'Not determined'", max_length=4000, null=True, blank=True)),
                ('published_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Type of end-point measurement: e.g. IC50, LD50, %inhibition etc, as it appears in the original publication', max_length=250, null=True, db_index=True, blank=True)),
                ('manual_curation_flag', chembl_core_db.db.customFields.ChemblPositiveIntegerField(default=0, null=True, choices=[(0, b'0'), (1, b'1'), (2, b'2')])),
                ('potential_duplicate', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Indicates whether the value is likely to be a repeat citation of a value reported in a previous ChEMBL paper, rather than a new, independent measurement.')),
                ('published_relation', chembl_core_db.db.customFields.ChemblCharField(help_text='Symbol constraining the activity value (e.g. >, <, =), as it appears in the original publication', max_length=50, null=True, db_index=True, blank=True)),
                ('original_activity_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(null=True, blank=True)),
                ('pchembl_value', models.DecimalField(decimal_places=2, max_digits=4, blank=True, help_text='Negative log of selected concentration-response activity values (IC50/EC50/XC50/AC50/Ki/Kd/Potency)', null=True, db_index=True)),
                ('bao_endpoint', chembl_core_db.db.customFields.ChemblCharField(help_text='ID for the corresponding result type in BioAssay Ontology (based on standard_type)', max_length=11, null=True, blank=True)),
                ('uo_units', chembl_core_db.db.customFields.ChemblCharField(help_text='ID for the corresponding unit in Unit Ontology (based on standard_units)', max_length=10, null=True, blank=True)),
                ('qudt_units', chembl_core_db.db.customFields.ChemblCharField(help_text='ID for the corresponding unit in QUDT Ontology (based on standard_units)', max_length=70, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table' : 'activities',   
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ActivityStdsLookup',
            fields=[
                ('std_act_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
                ('standard_type', chembl_core_db.db.customFields.ChemblCharField(help_text='The standard_type that other published_types in the activities table have been converted to.', max_length=250)),
                ('definition', chembl_core_db.db.customFields.ChemblCharField(help_text='A description/definition of the standard_type.', max_length=500, null=True, blank=True)),
                ('standard_units', chembl_core_db.db.customFields.ChemblCharField(help_text='The units that are applied to this standard_type and to which other published_units are converted. Note a standard_type may have more than one allowable standard_unit and therefore multiple rows in this table.', max_length=100)),
                ('normal_range_min', models.DecimalField(help_text="The lowest value for this activity type that is likely to be genuine. This is only an approximation, so lower genuine values may exist, but it may be desirable to validate these before using them. For a given standard_type/units, values in the activities table below this threshold are flagged with a data_validity_comment of 'Outside typical range'.", null=True, max_digits=24, decimal_places=12, blank=True)),
                ('normal_range_max', models.DecimalField(help_text="The highest value for this activity type that is likely to be genuine. This is only an approximation, so higher genuine values may exist, but it may be desirable to validate these before using them. For a given standard_type/units, values in the activities table above this threshold are flagged with a data_validity_comment of 'Outside typical range'.", null=True, max_digits=24, decimal_places=12, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table' : 'activity_stds_lookup'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='AssayParameters',
            fields=[
                ('assay_param_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Numeric primary key', serialize=False, primary_key=True)),
                ('parameter_value', chembl_core_db.db.customFields.ChemblCharField(help_text='The value of the particular parameter', max_length=2000)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'assay_parameters'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Assays',
            fields=[
                ('assay_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Unique ID for the assay', serialize=False, primary_key=True, blank=True)),
                ('description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of the reported assay', max_length=4000, null=True, db_index=True, blank=True)),
                ('assay_test_type', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=20, null=True, help_text='Type of assay system (i.e., in vivo or in vitro)', choices=[(b'In vivo', b'In vivo'), (b'In vitro', b'In vitro'), (b'Ex vivo', b'Ex vivo')])),
                ('assay_category', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=20, null=True, help_text='screening, confirmatory (ie: dose-response), summary, panel or other.', choices=[(b'screening', b'screening'), (b'panel', b'panel'), (b'confirmatory', b'confirmatory'), (b'summary', b'summary'), (b'other', b'other')])),
                ('assay_organism', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of the organism for the assay system (e.g., the organism, tissue or cell line in which an assay was performed). May differ from the target organism (e.g., for a human protein expressed in non-human cells, or pathogen-infected human cells).', max_length=250, null=True, blank=True)),
                ('assay_tax_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='NCBI tax ID for the assay organism.', null=True, blank=True)),
                ('assay_strain', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of specific strain of the assay organism used (where known)', max_length=200, null=True, blank=True)),
                ('assay_tissue', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of tissue used in the assay system (e.g., for tissue-based assays) or from which the assay system was derived (e.g., for cell/subcellular fraction-based assays).', max_length=100, null=True, blank=True)),
                ('assay_cell_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of cell type or cell line used in the assay system (e.g., for cell-based assays).', max_length=100, null=True, blank=True)),
                ('assay_subcellular_fraction', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of subcellular fraction used in the assay system (e.g., microsomes, mitochondria).', max_length=100, null=True, blank=True)),
                ('activity_count', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Number of activities recorded for this assay', null=True, blank=True)),
                ('assay_source', chembl_core_db.db.customFields.ChemblCharField(db_index=True, max_length=50, null=True, blank=True)),
                ('src_assay_id', chembl_core_db.db.customFields.ChemblCharField(help_text='Identifier for the assay in the source database/deposition (e.g., pubchem AID)', max_length=50, null=True, blank=True)),
                ('updated_on', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
                ('updated_by', chembl_core_db.db.customFields.ChemblCharField(max_length=250, null=True, blank=True)),
                ('orig_description', chembl_core_db.db.customFields.ChemblCharField(max_length=4000, null=True, blank=True)),
                ('a2t_complex', chembl_core_db.db.customFields.ChemblNullableBooleanField()),
                ('a2t_multi', chembl_core_db.db.customFields.ChemblNullableBooleanField()),
                ('mc_tax_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(null=True, blank=True)),
                ('mc_organism', chembl_core_db.db.customFields.ChemblCharField(max_length=100, null=True, blank=True)),
                ('mc_target_type', chembl_core_db.db.customFields.ChemblCharField(max_length=25, null=True, blank=True)),
                ('mc_target_name', chembl_core_db.db.customFields.ChemblCharField(max_length=4000, null=True, blank=True)),
                ('mc_target_accession', chembl_core_db.db.customFields.ChemblCharField(max_length=255, null=True, blank=True)),
                ('a2t_assay_tax_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(null=True, blank=True)),
                ('a2t_assay_organism', chembl_core_db.db.customFields.ChemblCharField(max_length=250, null=True, blank=True)),
                ('a2t_updated_on', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
                ('a2t_updated_by', chembl_core_db.db.customFields.ChemblCharField(max_length=100, null=True, blank=True)),
                ('bao_format', chembl_core_db.db.customFields.ChemblCharField(help_text='ID for the corresponding format type in BioAssay Ontology (e.g., cell-based, biochemical, organism-based etc)', max_length=11, null=True, db_index=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'assays'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='AssayType',
            fields=[
                ('assay_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Single character representing assay type', max_length=1, serialize=False, primary_key=True)),
                ('assay_desc', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of assay type', max_length=250, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'assay_type'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='AtcClassification',
            fields=[
                ('who_name', chembl_core_db.db.customFields.ChemblCharField(help_text='WHO/INN name for the compound', max_length=150, null=True, blank=True)),
                ('level1', chembl_core_db.db.customFields.ChemblCharField(help_text='First level of classification', max_length=10, null=True, blank=True)),
                ('level2', chembl_core_db.db.customFields.ChemblCharField(help_text='Second level of classification', max_length=10, null=True, blank=True)),
                ('level3', chembl_core_db.db.customFields.ChemblCharField(help_text='Third level of classification', max_length=10, null=True, blank=True)),
                ('level4', chembl_core_db.db.customFields.ChemblCharField(help_text='Fourth level of classification', max_length=10, null=True, blank=True)),
                ('level5', chembl_core_db.db.customFields.ChemblCharField(help_text='Complete ATC code for compound', max_length=10, serialize=False, primary_key=True)),
                ('who_id', chembl_core_db.db.customFields.ChemblCharField(help_text='WHO Identifier for compound', max_length=15, null=True, blank=True)),
                ('level1_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of first level of classification', max_length=150, null=True, blank=True)),
                ('level2_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of second level of classification', max_length=150, null=True, blank=True)),
                ('level3_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of third level of classification', max_length=150, null=True, blank=True)),
                ('level4_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of fourth level of classification', max_length=150, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table' : 'atc_classification'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='BindingSites',
            fields=[
                ('site_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key. Unique identifier for a binding site in a given target.', serialize=False, primary_key=True, blank=True)),
                ('site_name', chembl_core_db.db.customFields.ChemblCharField(help_text='Name/label for the binding site.', max_length=200, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'binding_sites'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='BioComponentSequences',
            fields=[
                ('component_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key. Unique identifier for each of the molecular components of biotherapeutics in ChEMBL (e.g., antibody chains, recombinant proteins, synthetic peptides).', serialize=False, primary_key=True, blank=True)),
                ('component_type', chembl_core_db.db.customFields.ChemblCharField(help_text="Type of molecular component (e.g., 'PROTEIN','DNA','RNA').", max_length=50)),
                ('description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description/name of molecular component.', max_length=200, null=True, blank=True)),
                ('sequence', chembl_core_db.db.customFields.ChemblTextField(help_text='Sequence of the biotherapeutic component.', null=True, blank=True)),
                ('sequence_md5sum', chembl_core_db.db.customFields.ChemblCharField(help_text='MD5 checksum of the sequence.', max_length=32, null=True, blank=True)),
                ('tax_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='NCBI tax ID for the species from which the sequence is derived. May be null for humanized monoclonal antibodies, synthetic peptides etc.', null=True, blank=True)),
                ('organism', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of the species from which the sequence is derived.', max_length=150, null=True, blank=True)),
                ('updated_on', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
                ('updated_by', chembl_core_db.db.customFields.ChemblCharField(max_length=100, null=True, blank=True)),
                ('insert_date', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
                ('accession', chembl_core_db.db.customFields.ChemblCharField(max_length=25, null=True, blank=True)),
                ('db_source', chembl_core_db.db.customFields.ChemblCharField(max_length=25, null=True, blank=True)),
                ('db_version', chembl_core_db.db.customFields.ChemblCharField(max_length=10, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'bio_component_sequences'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='BiotherapeuticComponents',
            fields=[
                ('biocomp_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'biotherapeutic_components'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CellDictionary',
            fields=[
                ('cell_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key. Unique identifier for each cell line in the target_dictionary.', serialize=False, primary_key=True, blank=True)),
                ('cell_name', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of each cell line (as used in the target_dicitonary pref_name).', max_length=50)),
                ('cell_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Longer description (where available) of the cell line.', max_length=200, null=True, blank=True)),
                ('cell_source_tissue', chembl_core_db.db.customFields.ChemblCharField(help_text='Tissue from which the cell line is derived, where known.', max_length=50, null=True, blank=True)),
                ('cell_source_organism', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of organism from which the cell line is derived.', max_length=150, null=True, blank=True)),
                ('cell_source_tax_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='NCBI tax ID of the organism from which the cell line is derived.', null=True, blank=True)),
                ('clo_id', chembl_core_db.db.customFields.ChemblCharField(help_text='ID for the corresponding cell line in Cell Line Ontology', max_length=11, null=True, blank=True)),
                ('efo_id', chembl_core_db.db.customFields.ChemblCharField(help_text='ID for the corresponding cell line in Experimental Factory Ontology', max_length=12, null=True, blank=True)),
                ('cellosaurus_id', chembl_core_db.db.customFields.ChemblCharField(help_text='ID for the corresponding cell line in Cellosaurus Ontology', max_length=15, null=True, blank=True)),
                ('downgraded', chembl_core_db.db.customFields.ChemblPositiveIntegerField(default=0, help_text='Indicates the cell line has been removed (if set to 1)', null=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'cell_dictionary'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ChemblIdLookup',
            fields=[
                ('chembl_id', chembl_core_db.db.customFields.ChemblCharField(help_text='ChEMBL identifier', max_length=20, serialize=False, primary_key=True)),
                ('entity_type', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=50, null=True, help_text='Type of entity (e.g., COMPOUND, ASSAY, TARGET)', choices=[(b'ASSAY', b'ASSAY'), (b'COMPOUND', b'COMPOUND'), (b'DOCUMENT', b'DOCUMENT'), (b'TARGET', b'TARGET')])),
                ('entity_id', chembl_core_db.db.customFields.ChemblIntegerField(help_text='Primary key for that entity in corresponding table (e.g., molregno for compounds, tid for targets)', null=True, blank=True)),
                ('status', chembl_core_db.db.customFields.ChemblCharField(default='ACTIVE', max_length=10, null=True, help_text='Indicates whether the status of the entity within the database - ACTIVE, INACTIVE (downgraded), OBS (obsolete/removed).', choices=[(b'ACTIVE', b'ACTIVE'), (b'INACTIVE', b'INACTIVE'), (b'OBS', b'OBS')])),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'chembl_id_lookup'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ComponentClass',
            fields=[
                ('comp_class_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'component_class'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ComponentDomains',
            fields=[
                ('compd_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
                ('start_position', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Start position of the domain within the sequence given in the component_sequences table.', null=True, blank=True)),
                ('end_position', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='End position of the domain within the sequence given in the component_sequences table.', null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'component_domains'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ComponentSequences',
            fields=[
                ('component_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key. Unique identifier for the component.', serialize=False, primary_key=True, blank=True)),
                ('component_type', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=50, null=True, help_text="Type of molecular component represented (e.g., 'PROTEIN','DNA','RNA').", choices=[(b'PROTEIN', b'PROTEIN'), (b'DNA', b'DNA'), (b'RNA', b'RNA')])),
                ('accession', chembl_core_db.db.customFields.ChemblCharField(help_text='Accession for the sequence in the source database from which it was taken (e.g., UniProt accession for proteins).', max_length=25, unique=True, null=True, blank=True)),
                ('sequence', chembl_core_db.db.customFields.ChemblTextField(help_text='A representative sequence for the molecular component, as given in the source sequence database (not necessarily the exact sequence used in the assay).', null=True, blank=True)),
                ('sequence_md5sum', chembl_core_db.db.customFields.ChemblCharField(help_text='MD5 checksum of the sequence.', max_length=32, null=True, blank=True)),
                ('description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description/name for the molecular component, usually taken from the source sequence database.', max_length=200, null=True, blank=True)),
                ('tax_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='NCBI tax ID for the sequence in the source database (i.e., species that the protein/nucleic acid sequence comes from).', null=True, blank=True)),
                ('organism', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of the organism the sequence comes from.', max_length=150, null=True, blank=True)),
                ('db_source', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=25, null=True, help_text='The name of the source sequence database from which sequences/accessions are taken. For UniProt proteins, this field indicates whether the sequence is from SWISS-PROT or TREMBL.', choices=[(b'Manual', b'Manual'), (b'SWISS-PROT', b'SWISS-PROT'), (b'TREMBL', b'TREMBL')])),
                ('db_version', chembl_core_db.db.customFields.ChemblCharField(help_text='The version of the source sequence database from which sequences/accession were last updated.', max_length=10, null=True, blank=True)),
                ('insert_date', chembl_core_db.db.customFields.ChemblDateField(default=datetime.date.today, null=True)),
                ('updated_on', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
                ('updated_by', chembl_core_db.db.customFields.ChemblCharField(max_length=100, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'component_sequences'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ComponentSynonyms',
            fields=[
                ('compsyn_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
                ('component_synonym', chembl_core_db.db.customFields.ChemblCharField(help_text='The synonym for the component.', max_length=500, null=True, blank=True)),
                ('syn_type', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=20, null=True, help_text='The type or origin of the synonym (e.g., GENE_SYMBOL).', choices=[(b'HGNC_SYMBOL', b'HGNC_SYMBOL'), (b'GENE_SYMBOL', b'GENE_SYMBOL'), (b'UNIPROT', b'UNIPROT'), (b'MANUAL', b'MANUAL'), (b'OTHER', b'OTHER'), (b'EC_NUMBER', b'EC_NUMBER')])),
                ('component', models.ForeignKey(help_text='Foreign key to the component_sequences table. The component to which this synonym applies.', to='chembl_core_model.ComponentSequences')),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'component_synonyms'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CompoundRecords',
            fields=[
                ('record_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Unique ID for a compound/record', serialize=False, primary_key=True, blank=True)),
                ('compound_key', chembl_core_db.db.customFields.ChemblCharField(help_text='Key text identifying this compound in the scientific document', max_length=250, null=True, db_index=True, blank=True)),
                ('compound_name', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of this compound recorded in the scientific document', max_length=4000, null=True, blank=True)),
                ('filename', chembl_core_db.db.customFields.ChemblCharField(max_length=250, null=True, blank=True)),
                ('old_compound_key', chembl_core_db.db.customFields.ChemblCharField(max_length=250, null=True, blank=True)),
                ('updated_by', chembl_core_db.db.customFields.ChemblCharField(max_length=100, null=True, blank=True)),
                ('updated_on', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
                ('src_compound_id', chembl_core_db.db.customFields.ChemblCharField(help_text='Identifier for the compound in the source database (e.g., pubchem SID)', max_length=150, null=True, db_index=True, blank=True)),
                ('removed', chembl_core_db.db.customFields.ChemblNullBooleanField(default=0)),
                ('src_compound_id_version', chembl_core_db.db.customFields.ChemblPositiveIntegerField(blank=True, null=True, choices=[(1, b'1'), (2, b'2'), (3, b'3'), (4, b'4'), (5, b'5'), (6, b'6'), (7, b'7'), (8, b'8')])),
                ('curated', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Can be marked as curated if the entry has been mapped to a molregno other than that given by the original structure, and hence care should be taken when updating')),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'compound_records'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ConfidenceScoreLookup',
            fields=[
                ('confidence_score', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='0-9 score showing level of confidence in assignment of the precise molecular target of the assay', serialize=False, primary_key=True)),
                ('description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of the target types assigned with each score', max_length=100)),
                ('target_mapping', chembl_core_db.db.customFields.ChemblCharField(help_text='Short description of the target types assigned with each score', max_length=30)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'confidence_score_lookup'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CurationLookup',
            fields=[
                ('curated_by', chembl_core_db.db.customFields.ChemblCharField(help_text='Short description of the level of curation', max_length=32, serialize=False, primary_key=True)),
                ('description', chembl_core_db.db.customFields.ChemblCharField(help_text='Definition of terms in the curated_by field.', max_length=100)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'curation_lookup'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='DataValidityLookup',
            fields=[
                ('data_validity_comment', chembl_core_db.db.customFields.ChemblCharField(help_text='Primary key. Short description of various types of errors/warnings applied to values in the activities table.', max_length=30, serialize=False, primary_key=True)),
                ('description', chembl_core_db.db.customFields.ChemblCharField(help_text='Definition of the terms in the data_validity_comment field.', max_length=200, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'data_validity_lookup'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='DefinedDailyDose',
            fields=[
                ('ddd_value', models.DecimalField(help_text='Value of defined daily dose', null=True, max_digits=9, decimal_places=2, blank=True)),
                ('ddd_units', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=20, null=True, help_text='Units of defined daily dose', choices=[(b'LSU', b'LSU'), (b'MU', b'MU'), (b'TU', b'TU'), (b'U', b'U'), (b'g', b'g'), (b'mcg', b'mcg'), (b'mg', b'mg'), (b'ml', b'ml'), (b'mmol', b'mmol'), (b'tablet', b'tablet')])),
                ('ddd_admr', chembl_core_db.db.customFields.ChemblCharField(help_text='Administration route for dose', max_length=30, null=True, blank=True)),
                ('ddd_comment', chembl_core_db.db.customFields.ChemblCharField(help_text='Comment', max_length=400, null=True, blank=True)),
                ('ddd_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Internal primary key', serialize=False, primary_key=True, blank=True)),
                ('atc_code', models.ForeignKey(db_column=b'atc_code', to='chembl_core_model.AtcClassification', help_text='ATC code for the compound (foreign key to ATC_CLASSIFICATION table)')),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'defined_daily_dose'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Docs',
            fields=[
                ('doc_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Unique ID for the document', serialize=False, primary_key=True, blank=True)),
                ('journal', chembl_core_db.db.customFields.ChemblCharField(help_text='Abbreviated journal name for an article', max_length=50, null=True, db_index=True, blank=True)),
                ('year', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Year of journal article publication', null=True, db_index=True, blank=True)),
                ('volume', chembl_core_db.db.customFields.ChemblCharField(help_text='Volume of journal article', max_length=50, null=True, db_index=True, blank=True)),
                ('issue', chembl_core_db.db.customFields.ChemblCharField(help_text='Issue of journal article', max_length=50, null=True, db_index=True, blank=True)),
                ('first_page', chembl_core_db.db.customFields.ChemblCharField(help_text='First page number of journal article', max_length=50, null=True, blank=True)),
                ('last_page', chembl_core_db.db.customFields.ChemblCharField(help_text='Last page number of journal article', max_length=50, null=True, blank=True)),
                ('pubmed_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='NIH pubmed record ID, where available', unique=True, null=True, blank=True)),
                ('updated_on', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
                ('updated_by', chembl_core_db.db.customFields.ChemblCharField(max_length=100, null=True, blank=True)),
                ('doi', chembl_core_db.db.customFields.ChemblCharField(help_text='Digital object identifier for this reference', max_length=50, null=True, blank=True)),
                ('title', chembl_core_db.db.customFields.ChemblCharField(help_text='Document title (e.g., Publication title or description of dataset)', max_length=500, null=True, blank=True)),
                ('doc_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Type of the document (e.g., Publication, Deposited dataset)', max_length=50, choices=[(b'PUBLICATION', b'PUBLICATION'), (b'BOOK', b'BOOK'), (b'DATASET', b'DATASET')])),
                ('authors', chembl_core_db.db.customFields.ChemblCharField(help_text='For a deposited dataset, the authors carrying out the screening and/or submitting the dataset.', max_length=4000, null=True, blank=True)),
                ('abstract', chembl_core_db.db.customFields.ChemblTextField(help_text='For a deposited dataset, a brief description of the dataset.', null=True, blank=True)),
                ('chembl', models.ForeignKey(blank=True, to='chembl_core_model.ChemblIdLookup', help_text='ChEMBL identifier for this document (for use on web interface etc)', unique=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'docs'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Domains',
            fields=[
                ('domain_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key. Unique identifier for each domain.', serialize=False, primary_key=True, blank=True)),
                ('domain_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Indicates the source of the domain (e.g., Pfam).', max_length=20, choices=[(b'Pfam-A', b'Pfam-A'), (b'Pfam-B', b'Pfam-B')])),
                ('source_domain_id', chembl_core_db.db.customFields.ChemblCharField(help_text='Identifier for the domain in the source database (e.g., Pfam ID such as PF00001).', max_length=20)),
                ('domain_name', chembl_core_db.db.customFields.ChemblCharField(help_text='Name given to the domain in the source database (e.g., 7tm_1).', max_length=20, null=True, blank=True)),
                ('domain_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Longer name or description for the domain.', max_length=500, null=True, blank=True)),
                ('component_sequences', models.ManyToManyField(to='chembl_core_model.ComponentSequences', null=True, through='chembl_core_model.ComponentDomains', blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'domains'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='DrugMechanism',
            fields=[
                ('mec_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key for each drug mechanism of action', serialize=False, primary_key=True, blank=True)),
                ('mechanism_of_action', chembl_core_db.db.customFields.ChemblCharField(help_text="Description of the mechanism of action e.g., 'Phosphodiesterase 5 inhibitor'", max_length=250, null=True, blank=True)),
                ('direct_interaction', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Flag to show whether the molecule is believed to interact directly with the target (1 = yes, 0 = no)')),
                ('molecular_mechanism', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Flag to show whether the mechanism of action describes the molecular target of the drug, rather than a higher-level physiological mechanism e.g., vasodilator (1 = yes, 0 = no)')),
                ('disease_efficacy', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Flag to show whether the target assigned is believed to play a role in the efficacy of the drug in the indication(s) for which it is approved (1 = yes, 0 = no)')),
                ('mechanism_comment', chembl_core_db.db.customFields.ChemblCharField(help_text='Additional comments regarding the mechanism of action', max_length=500, null=True, blank=True)),
                ('selectivity_comment', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=100, null=True, help_text='Additional comments regarding the selectivity of the drug', choices=[(b'Broad spectrum', b'Broad spectrum'), (b'EDG5 less relevant', b'EDG5 less relevant'), (b'M3 selective', b'M3 selective'), (b"Non-selective but type 5 receptor is overexpressed in Cushing's disease", b"Non-selective but type 5 receptor is overexpressed in Cushing's disease"), (b'Selective', b'Selective'), (b'Selective for the brain omega-1 receptor (i.e. BZ1-type, i.e. alpha1/beta1/gamma2-GABA receptor)', b'Selective for the brain omega-1 receptor (i.e. BZ1-type, i.e. alpha1/beta1/gamma2-GABA receptor)'), (b'Selectivity for types 2, 3 and 5', b'Selectivity for types 2, 3 and 5'), (b'selectivity for beta-3 containing complexes', b'selectivity for beta-3 containing complexes')])),
                ('binding_site_comment', chembl_core_db.db.customFields.ChemblCharField(help_text='Additional comments regarding the binding site of the drug', max_length=100, null=True, blank=True)),
                ('curated_by', chembl_core_db.db.customFields.ChemblCharField(max_length=20, null=True, blank=True)),
                ('date_added', chembl_core_db.db.customFields.ChemblDateField(default=datetime.date.today)),
                ('date_removed', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
                ('downgraded', chembl_core_db.db.customFields.ChemblNullableBooleanField()),
                ('downgrade_reason', chembl_core_db.db.customFields.ChemblCharField(max_length=200, null=True, blank=True)),
                ('uniprot_accessions', chembl_core_db.db.customFields.ChemblCharField(max_length=500, null=True, blank=True)),
                ('curator_comment', chembl_core_db.db.customFields.ChemblCharField(max_length=500, null=True, blank=True)),
                ('curation_status', chembl_core_db.db.customFields.ChemblCharField(default='PARTIAL', help_text='Show whether the curation for this row is complete', max_length=10, choices=[(b'COMPLETE', b'COMPLETE'), (b'PARTIAL', b'PARTIAL')])),
                ('action_type', models.ForeignKey(db_column=b'action_type', blank=True, to='chembl_core_model.ActionType', help_text='Type of action of the drug on the target e.g., agonist/antagonist etc (foreign key to action_type table)', null=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'drug_mechanism'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Formulations',
            fields=[
                ('ingredient', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of the approved ingredient within the product', max_length=200, null=True, blank=True)),
                ('strength', chembl_core_db.db.customFields.ChemblCharField(help_text='Dose strength', max_length=200, null=True, blank=True)),
                ('formulation_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
		'db_table': 'formulations'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='JournalArticles',
            fields=[
                ('int_pk', chembl_core_db.db.customFields.ChemblAutoField(serialize=False, primary_key=True, blank=True)),
                ('volume', chembl_core_db.db.customFields.ChemblPositiveIntegerField(db_index=True, null=True, blank=True)),
                ('issue', chembl_core_db.db.customFields.ChemblPositiveIntegerField(db_index=True, null=True, blank=True)),
                ('year', chembl_core_db.db.customFields.ChemblPositiveIntegerField(db_index=True, null=True, blank=True)),
                ('month', chembl_core_db.db.customFields.ChemblPositiveIntegerField(db_index=True, null=True, blank=True)),
                ('day', chembl_core_db.db.customFields.ChemblPositiveIntegerField(db_index=True, null=True, blank=True)),
                ('pagination', chembl_core_db.db.customFields.ChemblCharField(db_index=True, max_length=50, null=True, blank=True)),
                ('first_page', chembl_core_db.db.customFields.ChemblCharField(db_index=True, max_length=50, null=True, blank=True)),
                ('last_page', chembl_core_db.db.customFields.ChemblCharField(db_index=True, max_length=50, null=True, blank=True)),
                ('pubmed_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(db_index=True, null=True, blank=True)),
                ('doi', chembl_core_db.db.customFields.ChemblCharField(max_length=50, null=True, blank=True)),
                ('title', chembl_core_db.db.customFields.ChemblTextField(null=True, blank=True)),
                ('abstract', chembl_core_db.db.customFields.ChemblTextField(null=True, blank=True)),
                ('authors', chembl_core_db.db.customFields.ChemblTextField(null=True, blank=True)),
                ('year_raw', chembl_core_db.db.customFields.ChemblCharField(max_length=50, null=True, blank=True)),
                ('month_raw', chembl_core_db.db.customFields.ChemblCharField(max_length=50, null=True, blank=True)),
                ('day_raw', chembl_core_db.db.customFields.ChemblCharField(max_length=50, null=True, blank=True)),
                ('volume_raw', chembl_core_db.db.customFields.ChemblCharField(max_length=50, null=True, blank=True)),
                ('issue_raw', chembl_core_db.db.customFields.ChemblCharField(max_length=50, null=True, blank=True)),
                ('date_loaded', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'journal_articles'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Journals',
            fields=[
                ('journal_id', chembl_core_db.db.customFields.ChemblAutoField(serialize=False, primary_key=True, blank=True)),
                ('title', chembl_core_db.db.customFields.ChemblCharField(max_length=100, null=True, blank=True)),
                ('iso_abbreviation', chembl_core_db.db.customFields.ChemblCharField(max_length=50, null=True, blank=True)),
                ('issn_print', chembl_core_db.db.customFields.ChemblCharField(max_length=20, null=True, blank=True)),
                ('issn_electronic', chembl_core_db.db.customFields.ChemblCharField(max_length=20, null=True, blank=True)),
                ('publication_start_year', chembl_core_db.db.customFields.ChemblPositiveIntegerField(null=True, blank=True)),
                ('nlm_id', chembl_core_db.db.customFields.ChemblCharField(max_length=15, null=True, blank=True)),
                ('doc_journal', chembl_core_db.db.customFields.ChemblCharField(max_length=50, null=True, blank=True)),
                ('core_journal_flag', chembl_core_db.db.customFields.ChemblNullableBooleanField()),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'journals'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='LigandEff',
            fields=[
                ('activity', models.OneToOneField(primary_key=True, serialize=False, to='chembl_core_model.Activities', help_text='Link key to activities table')),
                ('bei', chembl_core_db.db.customFields.ChemblPositiveDecimalField(help_text='Binding Efficiency Index = p(XC50) *1000/MW_freebase', null=True, max_digits=9, decimal_places=2, blank=True)),
                ('sei', chembl_core_db.db.customFields.ChemblPositiveDecimalField(help_text='Surface Efficiency Index = p(XC50)*100/PSA', null=True, max_digits=9, decimal_places=2, blank=True)),
                ('le', chembl_core_db.db.customFields.ChemblPositiveDecimalField(help_text='Ligand Efficiency = deltaG/heavy_atoms  [from the Hopkins DDT paper 2004]', null=True, max_digits=9, decimal_places=2, blank=True)),
                ('lle', models.DecimalField(help_text='Lipophilic Ligand Efficiency = -logKi-ALogP. [from Leeson NRDD 2007]', null=True, max_digits=9, decimal_places=2, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'ligand_eff'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='MechanismRefs',
            fields=[
                ('mecref_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Primary key', serialize=False, primary_key=True)),
                ('ref_type', chembl_core_db.db.customFields.ChemblCharField(help_text="Type/source of reference (e.g., 'PubMed','DailyMed')", max_length=50, choices=[(b'ISBN', b'ISBN'), (b'IUPHAR', b'IUPHAR'), (b'DOI', b'DOI'), (b'EMA', b'EMA'), (b'PubMed', b'PUBMED'), (b'USPO', b'USPO'), (b'DailyMed', b'DAILYMED'), (b'FDA', b'FDA'), (b'Expert', b'EXPERT'), (b'Other', b'OTHER'), (b'InterPro', b'INTERPRO'), (b'Wikipedia', b'WIKIPEDIA'), (b'UniProt', b'UNIPROT'), (b'KEGG', b'KEGG')])),
                ('ref_id', chembl_core_db.db.customFields.ChemblCharField(help_text='Identifier for the reference in the source (e.g., PubMed ID or DailyMed setid)', max_length=100, null=True, blank=True)),
                ('ref_url', chembl_core_db.db.customFields.ChemblCharField(help_text='Full URL linking to the reference', max_length=200, null=True, blank=True)),
                ('mechanism', models.ForeignKey(db_column=b'mec_id', to='chembl_core_model.DrugMechanism', help_text='Foreign key to drug_mechanism table - indicating the mechanism to which the references refer')),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'mechanism_refs'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='MoleculeAtcClassification',
            fields=[
                ('mol_atc_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key', serialize=False, primary_key=True, blank=True)),
                ('atc_classification', models.ForeignKey(db_column=b'level5', to='chembl_core_model.AtcClassification', help_text='ATC code (foreign key to atc_classification table)')),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'molecule_atc_classification'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='MoleculeDictionary',
            fields=[
                ('molregno', chembl_core_db.db.customFields.ChemblAutoField(help_text='Internal Primary Key for the molecule', serialize=False, primary_key=True, blank=True)),
                ('pref_name', chembl_core_db.db.customFields.ChemblCharField(help_text='Preferred name for the molecule', max_length=255, null=True, db_index=True, blank=True)),
                ('max_phase', chembl_core_db.db.customFields.ChemblPositiveIntegerField(default=0, help_text='Maximum phase of development reached for the compound (4 = approved). Null where max phase has not yet been assigned.', db_index=True, choices=[(0, b'0'), (1, b'1'), (2, b'2'), (3, b'3'), (4, b'4')])),
                ('therapeutic_flag', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Indicates that a drug has a therapeutic application (as opposed to e.g., an imaging agent, additive etc).', db_index=True)),
                ('dosed_ingredient', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Indicates that the drug is dosed in this form (e.g., a particular salt)')),
                ('structure_key', chembl_core_db.db.customFields.ChemblCharField(null=True, max_length=27, blank=True, help_text='Unique key for the structure/sequence (e.g., inchi_key or sequence md5sum) to help enforce non-redundancy.', unique=True, db_index=True)),
                ('structure_type', chembl_core_db.db.customFields.ChemblCharField(default=b'MOL', help_text='Indications whether the molecule has a small molecule structure or a protein sequence (MOL indicates an entry in the compound_structures table, SEQ indications an entry in the protein_therapeutics table, NONE indicates an entry in neither table, e.g., structure unknown)', max_length=10, choices=[(b'NONE', b'NONE'), (b'MOL', b'MOL'), (b'SEQ', b'SEQ'), (b'BOTH', b'BOTH')])),
                ('chebi_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Assigned ChEBI ID for the compound, where it is a small molecule.', unique=True, null=True, blank=True)),
                ('chebi_par_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Preferred ChEBI ID for the compound (where different from assigned)', null=True, blank=True)),
                ('insert_date', chembl_core_db.db.customFields.ChemblDateField(default=datetime.date.today, null=True)),
                ('molfile_update', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
                ('downgraded', chembl_core_db.db.customFields.ChemblBooleanField(default=False)),
                ('downgrade_reason', chembl_core_db.db.customFields.ChemblCharField(max_length=2000, null=True, blank=True)),
                ('replacement_mrn', chembl_core_db.db.customFields.ChemblPositiveIntegerField(null=True, blank=True)),
                ('checked_by', chembl_core_db.db.customFields.ChemblCharField(max_length=2000, null=True, blank=True)),
                ('nomerge', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text="Flag to show that this entry shouldn't be merged with others of the same structure (when set to 1)")),
                ('nomerge_reason', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=200, null=True, help_text='Reason for entry not being merged with others of the same structure (e.g., known to be a stereoisomer)', choices=[(b'GSK', b'GSK'), (b'PARENT', b'PARENT'), (b'PDBE', b'PDBE'), (b'SALT', b'SALT')])),
                ('molecule_type', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=30, null=True, help_text='Type of molecule (Small molecule, Protein, Antibody, Oligosaccharide, Oligonucleotide, Cell, Unknown)', choices=[(b'Antibody', b'Antibody'), (b'Cell', b'Cell'), (b'Enzyme', b'Enzyme'), (b'Oligonucleotide', b'Oligonucleotide'), (b'Oligosaccharide', b'Oligosaccharide'), (b'Protein', b'Protein'), (b'Small molecule', b'Small molecule'), (b'Unclassified', b'Unclassified'), (b'Unknown', b'Unknown')])),
                ('first_approval', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Earliest known approval year for the molecule', null=True, blank=True)),
                ('oral', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Indicates whether the drug is known to be administered orally.')),
                ('parenteral', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Indicates whether the drug is known to be administered parenterally')),
                ('topical', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Indicates whether the drug is known to be administered topically.')),
                ('black_box_warning', chembl_core_db.db.customFields.ChemblNullBooleanField(default=0, help_text='Indicates that the drug has a black box warning')),
                ('natural_product', chembl_core_db.db.customFields.ChemblNullBooleanField(default=-1, help_text='Indicates whether the compound is natural product-derived (currently curated only for drugs)')),
                ('first_in_class', chembl_core_db.db.customFields.ChemblNullBooleanField(default=-1, help_text='Indicates whether this is known to be the first compound of its class (e.g., acting on a particular target).')),
                ('chirality', chembl_core_db.db.customFields.ChemblIntegerField(default=-1, help_text='Shows whether a drug is dosed as a racemic mixture (0), single stereoisomer (1) or is an achiral molecule (2)', choices=[(-1, b'-1'), (0, b'0'), (1, b'1'), (2, b'2')])),
                ('prodrug', chembl_core_db.db.customFields.ChemblNullBooleanField(default=-1, help_text='Indicates that the molecule is a pro-drug (see molecule hierarchy for active component, where known)')),
                ('exclude', chembl_core_db.db.customFields.ChemblBooleanField(default=False)),
                ('inorganic_flag', chembl_core_db.db.customFields.ChemblNullBooleanField(default=0, help_text='Indicates whether the molecule is inorganic (i.e., containing only metal atoms and <2 carbon atoms)')),
                ('usan_year', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='The year in which the application for a USAN/INN name was made', null=True, blank=True)),
                ('availability_type', chembl_core_db.db.customFields.ChemblIntegerField(blank=True, help_text='The availability type for the drug (0 = discontinued, 1 = prescription only, 2 = over the counter)', null=True, choices=[(-1, b'-1'), (0, b'0'), (1, b'1'), (2, b'2')])),
                ('usan_stem', chembl_core_db.db.customFields.ChemblCharField(help_text='Where the compound has been assigned a USAN name, this indicates the stem, as described in the USAN_STEM table.', max_length=50, null=True, blank=True)),
                ('polymer_flag', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Indicates whether a molecule is a small molecule polymer (e.g., polistyrex)')),
                ('usan_substem', chembl_core_db.db.customFields.ChemblCharField(help_text='Where the compound has been assigned a USAN name, this indicates the substem', max_length=50, null=True, blank=True)),
                ('usan_stem_definition', chembl_core_db.db.customFields.ChemblCharField(help_text='Definition of the USAN stem', max_length=1000, null=True, blank=True)),
                ('indication_class', chembl_core_db.db.customFields.ChemblCharField(help_text='Indication class(es) assigned to a drug in the USP dictionary', max_length=1000, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'molecule_dictionary'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CompoundStructures',
            fields=[
                ('molecule', models.OneToOneField(primary_key=True, db_column=b'molregno', serialize=False, to='chembl_core_model.MoleculeDictionary', help_text='Internal Primary Key for the compound structure and foreign key to molecule_dictionary table')),
                ('molfile', chembl_core_db.db.customFields.ChemblTextField(help_text='MDL Connection table representation of compound', null=True, blank=True)),
                ('standard_inchi', chembl_core_db.db.customFields.ChemblCharField(null=True, max_length=4000, blank=True, help_text='IUPAC standard InChI for the compound', unique=True, db_index=True)),
                ('standard_inchi_key', chembl_core_db.db.customFields.ChemblCharField(help_text='IUPAC standard InChI key for the compound', max_length=27, db_index=True)),
                ('canonical_smiles', chembl_core_db.db.customFields.ChemblCharField(help_text='Canonical smiles, generated using pipeline pilot', max_length=4000, null=True, db_index=True, blank=True)),
                ('structure_exclude_flag', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Indicates whether the structure for this compound should be hidden from users (e.g., organometallic compounds with bad valence etc)')),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'compound_structures'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CompoundProperties',
            fields=[
                ('molecule', models.OneToOneField(primary_key=True, db_column=b'molregno', serialize=False, to='chembl_core_model.MoleculeDictionary', help_text='Foreign key to compounds table (compound structure)')),
                ('mw_freebase', chembl_core_db.db.customFields.ChemblPositiveDecimalField(decimal_places=2, max_digits=9, blank=True, help_text='Molecular weight of parent compound', null=True, db_index=True)),
                ('alogp', models.DecimalField(decimal_places=2, max_digits=9, blank=True, help_text='Calculated ALogP', null=True, db_index=True)),
                ('hba', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Number hydrogen bond acceptors', null=True, db_index=True, blank=True)),
                ('hbd', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Number hydrogen bond donors', null=True, db_index=True, blank=True)),
                ('psa', chembl_core_db.db.customFields.ChemblPositiveDecimalField(decimal_places=2, max_digits=9, blank=True, help_text='Polar surface area', null=True, db_index=True)),
                ('rtb', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Number rotatable bonds', null=True, db_index=True, blank=True)),
                ('ro3_pass', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=3, null=True, help_text='Indicates whether the compound passes the rule-of-three (mw < 300, logP < 3 etc)', choices=[(b'Y', b'Yes'), (b'N', b'No')])),
                ('num_ro5_violations', chembl_core_db.db.customFields.ChemblPositiveIntegerField(blank=True, help_text='Number of violations of rule-of-five', null=True, db_index=True, choices=[(0, b'0'), (1, b'1'), (2, b'2'), (3, b'3'), (4, b'4')])),
                ('med_chem_friendly', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=3, null=True, help_text='Indicates whether the compound is considered Med Chem friendly (Y/N)', choices=[(b'Y', b'Yes'), (b'N', b'No')])),
                ('acd_most_apka', models.DecimalField(help_text='The most acidic pKa calculated using ACDlabs v12.01', null=True, max_digits=9, decimal_places=2, blank=True)),
                ('acd_most_bpka', models.DecimalField(help_text='The most basic pKa calculated using ACDlabs v12.01', null=True, max_digits=9, decimal_places=2, blank=True)),
                ('acd_logp', models.DecimalField(help_text='The calculated octanol/water partition coefficient using ACDlabs v12.01', null=True, max_digits=9, decimal_places=2, blank=True)),
                ('acd_logd', models.DecimalField(help_text='The calculated octanol/water distribution coefficient at pH7.4 using ACDlabs v12.01', null=True, max_digits=9, decimal_places=2, blank=True)),
                ('molecular_species', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=50, null=True, help_text='Indicates whether the compound is an acid/base/neutral', choices=[(b'ACID', b'ACID'), (b'BASE', b'BASE'), (b'ZWITTERION', b'ZWITTERION'), (b'NEUTRAL', b'NEUTRAL')])),
                ('full_mwt', chembl_core_db.db.customFields.ChemblPositiveDecimalField(help_text='Molecular weight of the full compound including any salts', null=True, max_digits=9, decimal_places=2, blank=True)),
                ('aromatic_rings', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Number of aromatic rings', null=True, blank=True)),
                ('heavy_atoms', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Number of heavy (non-hydrogen) atoms', null=True, blank=True)),
                ('num_alerts', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Number of structural alerts (as defined by Brenk et al., ChemMedChem 2008)', null=True, blank=True)),
                ('qed_weighted', chembl_core_db.db.customFields.ChemblPositiveDecimalField(help_text='Weighted quantitative estimate of drug likeness (as defined by Bickerton et al., Nature Chem 2012)', null=True, max_digits=3, decimal_places=2, blank=True)),
                ('updated_on', chembl_core_db.db.customFields.ChemblDateField(help_text='Shows date properties were last recalculated', null=True, blank=True)),
                ('mw_monoisotopic', chembl_core_db.db.customFields.ChemblPositiveDecimalField(help_text='Monoisotopic parent molecular weight', null=True, max_digits=11, decimal_places=4, blank=True)),
                ('full_molformula', chembl_core_db.db.customFields.ChemblCharField(help_text='Molecular formula for the full compound (including any salt)', max_length=100, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'compound_properties'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CompoundMols',
            fields=[
                ('molecule', models.OneToOneField(primary_key=True, db_column=b'molregno', serialize=False, to='chembl_core_model.MoleculeDictionary')),
                ('ctab', chembl_core_db.db.customFields.BlobField(null=True, db_column=b'ctab', blank=True)),
            ],
            options={
                'abstract': False,
                'db_table': 'compound_mols',
                'managed': True,
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CompoundImages',
            fields=[
                ('molecule', models.OneToOneField(primary_key=True, db_column=b'molregno', serialize=False, to='chembl_core_model.MoleculeDictionary')),
                ('png', chembl_core_db.db.customFields.BlobField(null=True, blank=True)),
                ('png_500', chembl_core_db.db.customFields.BlobField(null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'compound_images'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Biotherapeutics',
            fields=[
                ('molecule', models.OneToOneField(primary_key=True, db_column=b'molregno', serialize=False, to='chembl_core_model.MoleculeDictionary', help_text='Foreign key to molecule_dictionary')),
                ('description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of the biotherapeutic.', max_length=2000, null=True, blank=True)),
                ('bio_component_sequences', models.ManyToManyField(to='chembl_core_model.BioComponentSequences', null=True, through='chembl_core_model.BiotherapeuticComponents', blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'biotherapeutics'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='MoleculeHierarchy',
            fields=[
                ('molecule', models.OneToOneField(primary_key=True, db_column=b'molregno', serialize=False, to='chembl_core_model.MoleculeDictionary', help_text='Foreign key to compounds table. This field holds a list of all of the ChEMBL compounds with associated data (e.g., activity information, approved drugs). Parent compounds that are generated only by removing salts, and which do not themselves have any associated data will not appear here.')),
                ('active_molecule', models.ForeignKey(related_name='active', db_column=b'active_molregno', blank=True, to='chembl_core_model.MoleculeDictionary', help_text="Where a compound is a pro-drug, this represents the active metabolite of the 'dosed' compound given by parent_molregno. Where parent_molregno and active_molregno are the same, the compound is not currently known to be a pro-drug. ", null=True)),
                ('parent_molecule', models.ForeignKey(related_name='parent', db_column=b'parent_molregno', blank=True, to='chembl_core_model.MoleculeDictionary', help_text='Represents parent compound of molregno in first field (i.e., generated by removing salts). Where molregno and parent_molregno are same, the initial ChEMBL compound did not contain a salt component, or else could not be further processed for various reasons (e.g., inorganic mixture). Compounds which are only generated by removing salts will appear in this field only. Those which, themselves, have any associated data (e.g., activity data) or are launched drugs will also appear in the molregno field.', null=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'molecule_hierarchy'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='MoleculeSynonyms',
            fields=[
                ('syn_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Type of name/synonym (e.g., TRADE_NAME, RESEARCH_CODE, USAN)', max_length=50)),
                ('molsyn_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
                ('synonyms', chembl_core_db.db.customFields.ChemblCharField(help_text='Synonym for the compound', max_length=200, null=True, blank=True)),
                ('molecule', models.ForeignKey(db_column=b'molregno', to='chembl_core_model.MoleculeDictionary', help_text='Foreign key to molecule_dictionary')),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'molecule_synonyms'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='OrganismClass',
            fields=[
                ('oc_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Internal primary key', serialize=False, primary_key=True, blank=True)),
                ('tax_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='NCBI taxonomy ID for the organism (corresponding to tax_ids in assay2target and target_dictionary tables)', unique=True, null=True, blank=True)),
                ('l1', chembl_core_db.db.customFields.ChemblCharField(help_text='Highest level classification (e.g., Eukaryotes, Bacteria, Fungi etc)', max_length=200, null=True, blank=True)),
                ('l2', chembl_core_db.db.customFields.ChemblCharField(help_text='Second level classification', max_length=200, null=True, blank=True)),
                ('l3', chembl_core_db.db.customFields.ChemblCharField(help_text='Third level classification', max_length=200, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'organism_class'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ParameterType',
            fields=[
                ('parameter_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Short name for the type of parameter associated with an assay', max_length=20, serialize=False, primary_key=True)),
                ('description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of the parameter type', max_length=2000, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'parameter_type'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PredictedBindingDomains',
            fields=[
                ('predbind_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
                ('prediction_method', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=50, null=True, help_text="The method used to assign the binding domain (e.g., 'Single domain' where the protein has only 1 domain, 'Multi domain' where the protein has multiple domains, but only 1 is known to bind small molecules in other proteins).", choices=[(b'Manual', b'Manual'), (b'Multi domain', b'Multi domain'), (b'Single domain', b'Single domain')])),
                ('confidence', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=10, null=True, help_text='The level of confidence assigned to the prediction (high where the protein has only 1 domain, medium where the compound has multiple domains, but only 1 known small molecule-binding domain).', choices=[(b'high', b'high'), (b'medium', b'medium'), (b'low', b'low')])),
                ('activity', models.ForeignKey(blank=True, to='chembl_core_model.Activities', help_text='Foreign key to the activities table, indicating the compound/assay(+target) combination for which this prediction is made.', null=True)),
                ('site', models.ForeignKey(blank=True, to='chembl_core_model.BindingSites', help_text='Foreign key to the binding_sites table, indicating the binding site (domain) that the compound is predicted to bind to.', null=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'predicted_binding_domains'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Products',
            fields=[
                ('dosage_form', chembl_core_db.db.customFields.ChemblCharField(help_text='The dosage form of the product (e.g., tablet, capsule etc)', max_length=200, null=True, blank=True)),
                ('route', chembl_core_db.db.customFields.ChemblCharField(help_text='The administration route of the product (e.g., oral, injection etc)', max_length=200, null=True, blank=True)),
                ('trade_name', chembl_core_db.db.customFields.ChemblCharField(help_text='The trade name for the product', max_length=200, null=True, blank=True)),
                ('approval_date', chembl_core_db.db.customFields.ChemblDateField(help_text='The FDA approval date for the product (not necessarily first approval of the active ingredient)', null=True, blank=True)),
                ('ad_type', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=5, null=True, help_text='RX = prescription, OTC = over the counter, DISCN = discontinued', choices=[(b'OTC', b'OTC'), (b'RX', b'RX'), (b'DISCN', b'DISCN')])),
                ('oral', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Flag to show whether product is orally delivered')),
                ('topical', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Flag to show whether product is topically delivered')),
                ('parenteral', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Flag to show whether product is parenterally delivered')),
                ('information_source', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=100, null=True, help_text='Source of the product information (e.g., Orange Book)', choices=[(b'CBER', b'CBER'), (b'CDER', b'CDER'), (b'MANUAL', b'MANUAL'), (b'ORANGE BOOK', b'ORANGE BOOK')])),
                ('black_box_warning', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Flag to show whether the product label has a black box warning')),
                ('product_class', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=30, null=True, choices=[(b'VACCINE', b'Vaccine'), (b'ANTI-RHESIS ANTIBODY', b'Anti-rhesis antibody')])),
                ('applicant_full_name', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of the company applying for FDA approval', max_length=200, null=True, blank=True)),
                ('innovator_company', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Flag to show whether the applicant is the innovator of the product')),
                ('product_id', chembl_core_db.db.customFields.ChemblCharField(help_text='FDA application number for the product', max_length=30, serialize=False, primary_key=True)),
                ('load_date', chembl_core_db.db.customFields.ChemblDateField(help_text='The date on which one or more of the following fields were created or updated: doasge_form, route, trade_name, approval_date, ad_type, oral, topical, parenteral, information_source, or applicant_full_name). This date is assigned by the EBI parser.', null=True, blank=True)),
                ('removed_date', chembl_core_db.db.customFields.ChemblDateField(help_text="The date on which this product was first identified (by ebi parser) as having been removed from the information source. The recording of this date was first initiated on 30-JUN-10. Note that a small number of products are removed from OB, but then re-appear... in these cases this field is re-set to 'null' when the product re-appears..", null=True, blank=True)),
                ('nda_type', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=10, null=True, help_text='New Drug Application Type. The type of new drug application approval.  New Drug Applications (NDA or innovator)  are "N".   Abbreviated New Drug Applications (ANDA or generic) are "A".', choices=[(b'A', b'A'), (b'N', b'N')])),
                ('tmp_ingred_count', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Number of ingredients in the product, to show which are combinations', null=True, blank=True)),
                ('exclude', chembl_core_db.db.customFields.ChemblIntegerField(help_text='Non-FDA products, to be excluded', null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True, 'db_table': 'products'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinClassification',
            fields=[
                ('protein_class_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key. Unique identifier for each protein family classification.', serialize=False, primary_key=True, blank=True)),
                ('parent_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Protein_class_id for the parent of this protein family.', null=True, blank=True)),
                ('pref_name', chembl_core_db.db.customFields.ChemblCharField(help_text='Preferred/full name for this protein family.', max_length=500, null=True, blank=True)),
                ('short_name', chembl_core_db.db.customFields.ChemblCharField(help_text='Short/abbreviated name for this protein family (not necessarily unique).', max_length=50, null=True, blank=True)),
                ('protein_class_desc', chembl_core_db.db.customFields.ChemblCharField(help_text='Concatenated description of each classification for searching purposes etc.', max_length=410)),
                ('definition', chembl_core_db.db.customFields.ChemblCharField(help_text='Definition of the protein family.', max_length=4000, null=True, blank=True)),
                ('downgraded', chembl_core_db.db.customFields.ChemblNullableBooleanField()),
                ('replaced_by', chembl_core_db.db.customFields.ChemblPositiveIntegerField(null=True, blank=True)),
                ('class_level', chembl_core_db.db.customFields.ChemblPositiveIntegerField(blank=True, null=True, choices=[(0, b'0'), (1, b'1'), (2, b'2'), (3, b'3'), (4, b'4'), (5, b'5'), (6, b'6'), (7, b'7'), (8, b'8')])),
                ('sort_order', chembl_core_db.db.customFields.ChemblPositiveIntegerField(null=True, blank=True)),
                ('component_sequences', models.ManyToManyField(to='chembl_core_model.ComponentSequences', null=True, through='chembl_core_model.ComponentClass', blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'protein_classification'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinClassSynonyms',
            fields=[
                ('protclasssyn_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
                ('protein_class_synonym', chembl_core_db.db.customFields.ChemblCharField(help_text='The synonym for the protein class.', max_length=1000, null=True, blank=True)),
                ('syn_type', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=20, null=True, help_text='The type or origin of the synonym (e.g., ChEMBL, Concept Wiki, UMLS).', choices=[(b'CHEMBL', b'CHEMBL'), (b'CONCEPT_WIKI', b'CONCEPT_WIKI'), (b'UMLS', b'UMLS'), (b'CW_XREF', b'CW_XREF'), (b'MESH_XREF', b'MESH_XREF')])),
                ('protein_class', models.ForeignKey(help_text='Foreign key to the PROTEIN_CLASSIFICATION table. The protein_class to which this synonym applies.', to='chembl_core_model.ProteinClassification')),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'protein_class_synonyms'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProteinFamilyClassification',
            fields=[
                ('protein_class_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key. Unique identifier for each classification.', serialize=False, primary_key=True, blank=True)),
                ('protein_class_desc', chembl_core_db.db.customFields.ChemblCharField(help_text='Concatenated description of each classification for searching purposes etc.', unique=True, max_length=810)),
                ('l1', chembl_core_db.db.customFields.ChemblCharField(help_text='First level classification (e.g., Enzyme, Transporter, Ion Channel).', max_length=100)),
                ('l2', chembl_core_db.db.customFields.ChemblCharField(help_text='Second level classification.', max_length=100, null=True, blank=True)),
                ('l3', chembl_core_db.db.customFields.ChemblCharField(help_text='Third level classification.', max_length=100, null=True, blank=True)),
                ('l4', chembl_core_db.db.customFields.ChemblCharField(help_text='Fourth level classification.', max_length=100, null=True, blank=True)),
                ('l5', chembl_core_db.db.customFields.ChemblCharField(help_text='Fifth level classification.', max_length=100, null=True, blank=True)),
                ('l6', chembl_core_db.db.customFields.ChemblCharField(help_text='Sixth level classification.', max_length=100, null=True, blank=True)),
                ('l7', chembl_core_db.db.customFields.ChemblCharField(help_text='Seventh level classification.', max_length=100, null=True, blank=True)),
                ('l8', chembl_core_db.db.customFields.ChemblCharField(help_text='Eighth level classification.', max_length=100, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'protein_family_classification'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='RecordDrugProperties',
            fields=[
                ('record', models.OneToOneField(primary_key=True, serialize=False, to='chembl_core_model.CompoundRecords')),
                ('max_phase', chembl_core_db.db.customFields.ChemblPositiveIntegerField(default=0, help_text='Maximum phase of development reached for the compound (4 = approved). Null where max phase has not yet been assigned.', db_index=True, choices=[(0, b'0'), (1, b'1'), (2, b'2'), (3, b'3'), (4, b'4')])),
                ('withdrawn_status', chembl_core_db.db.customFields.ChemblCharField(max_length=10, null=True, blank=True)),
                ('molecule_type', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=30, null=True, help_text='Type of molecule (Small molecule, Protein, Antibody, Oligosaccharide, Oligonucleotide, Cell, Unknown)', choices=[(b'Antibody', b'Antibody'), (b'Cell', b'Cell'), (b'Enzyme', b'Enzyme'), (b'Oligonucleotide', b'Oligonucleotide'), (b'Oligosaccharide', b'Oligosaccharide'), (b'Protein', b'Protein'), (b'Small molecule', b'Small molecule'), (b'Unclassified', b'Unclassified'), (b'Unknown', b'Unknown')])),
                ('first_approval', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Earliest known approval year for the molecule', null=True, blank=True)),
                ('oral', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Indicates whether the drug is known to be administered orally.')),
                ('parenteral', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Indicates whether the drug is known to be administered parenterally')),
                ('topical', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Indicates whether the drug is known to be administered topically.')),
                ('black_box_warning', chembl_core_db.db.customFields.ChemblNullBooleanField(default=0, help_text='Indicates that the drug has a black box warning')),
                ('first_in_class', chembl_core_db.db.customFields.ChemblNullBooleanField(default=-1, help_text='Indicates whether this is known to be the first compound of its class (e.g., acting on a particular target).')),
                ('chirality', chembl_core_db.db.customFields.ChemblIntegerField(default=-1, help_text='Shows whether a drug is dosed as a racemic mixture (0), single stereoisomer (1) or is an achiral molecule (2)', choices=[(-1, b'-1'), (0, b'0'), (1, b'1'), (2, b'2')])),
                ('prodrug', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Indicates that the molecule is a pro-drug (see molecule hierarchy for active component, where known)')),
                ('therapeutic_flag', chembl_core_db.db.customFields.ChemblBooleanField(default=False, help_text='Indicates that a drug has a therapeutic application (as opposed to e.g., an imaging agent, additive etc).', db_index=True)),
                ('natural_product', chembl_core_db.db.customFields.ChemblNullBooleanField(default=-1, help_text='Indicates whether the compound is natural product-derived (currently curated only for drugs)')),
                ('inorganic_flag', chembl_core_db.db.customFields.ChemblNullBooleanField(default=0, help_text='Indicates whether the molecule is inorganic (i.e., containing only metal atoms and <2 carbon atoms)')),
                ('applicants', chembl_core_db.db.customFields.ChemblCharField(max_length=1000, null=True, blank=True)),
                ('usan_stem', chembl_core_db.db.customFields.ChemblCharField(help_text='Where the compound has been assigned a USAN name, this indicates the stem, as described in the USAN_STEM table.', max_length=50, null=True, blank=True)),
                ('usan_year', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='The year in which the application for a USAN/INN name was made', null=True, blank=True)),
                ('availability_type', chembl_core_db.db.customFields.ChemblIntegerField(blank=True, help_text='The availability type for the drug (0 = discontinued, 1 = prescription only, 2 = over the counter)', null=True, choices=[(-1, b'-1'), (0, b'0'), (1, b'1'), (2, b'2')])),
                ('usan_substem', chembl_core_db.db.customFields.ChemblCharField(help_text='Where the compound has been assigned a USAN name, this indicates the substem', max_length=50, null=True, blank=True)),
                ('indication_class', chembl_core_db.db.customFields.ChemblCharField(help_text='Indication class(es) assigned to a drug in the USP dictionary', max_length=1000, null=True, blank=True)),
                ('usan_stem_definition', chembl_core_db.db.customFields.ChemblCharField(help_text='Definition of the USAN stem', max_length=1000, null=True, blank=True)),
                ('polymer_flag', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text='Indicates whether a molecule is a small molecule polymer (e.g., polistyrex)')),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table': 'record_drug_properties'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='RelationshipType',
            fields=[
                ('relationship_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Relationship_type flag used in the assay2target table', max_length=1, serialize=False, primary_key=True)),
                ('relationship_desc', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of relationship_type flags', max_length=250, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
                'db_table' : 'relationship_type'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ResearchCompanies',
            fields=[
                ('co_stem_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
                ('company', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of current company associated with this research code stem.', max_length=100, null=True, blank=True)),
                ('country', chembl_core_db.db.customFields.ChemblCharField(help_text='Country in which the company uses this research code stem.', max_length=50, null=True, blank=True)),
                ('previous_company', chembl_core_db.db.customFields.ChemblCharField(help_text='Previous name of the company associated with this research code stem (e.g., if the company has undergone acquisitions/mergers).', max_length=100, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
'db_table': 'research_companies'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ResearchStem',
            fields=[
                ('res_stem_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key. Unique ID for each research code stem.', serialize=False, primary_key=True, blank=True)),
                ('research_stem', chembl_core_db.db.customFields.ChemblCharField(help_text='The actual stem/prefix used in the research code.', max_length=20, unique=True, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
            'db_table': 'research_stem'},
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='SiteComponents',
            fields=[
                ('sitecomp_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
                ('site_residues', chembl_core_db.db.customFields.ChemblCharField(help_text='List of residues from the given molecular component that make up the binding site (where not know, will be null).', max_length=2000, null=True, blank=True)),
                ('component', models.ForeignKey(blank=True, to='chembl_core_model.ComponentSequences', help_text='Foreign key to the component_sequences table, indicating which molecular component of the target is involved in the binding site.', null=True)),
                ('domain', models.ForeignKey(blank=True, to='chembl_core_model.Domains', help_text='Foreign key to the domains table, indicating which domain of the given molecular component is involved in the binding site (where not known, the domain_id may be null).', null=True)),
                ('site', models.ForeignKey(help_text='Foreign key to binding_sites table.', to='chembl_core_model.BindingSites')),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'site_components'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Source',
            fields=[
                ('src_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Identifier for each source (used in compound_records and assays tables)', serialize=False, primary_key=True, blank=True)),
                ('src_description', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of the data source', max_length=500, null=True, blank=True)),
                ('src_short_name', chembl_core_db.db.customFields.ChemblCharField(help_text='A short name for each data source, for display purposes', max_length=20, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'source'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TargetComponents',
            fields=[
                ('relationship', chembl_core_db.db.customFields.ChemblCharField(default='SUBUNIT', max_length=20, choices=[(b'COMPARATIVE PROTEIN', b'COMPARATIVE PROTEIN'), (b'EQUIVALENT PROTEIN', b'EQUIVALENT PROTEIN'), (b'FUSION PROTEIN', b'FUSION PROTEIN'), (b'GROUP MEMBER', b'GROUP MEMBER'), (b'INTERACTING PROTEIN', b'INTERACTING PROTEIN'), (b'PROTEIN SUBUNIT', b'PROTEIN SUBUNIT'), (b'RNA', b'RNA'), (b'RNA SUBUNIT', b'RNA SUBUNIT'), (b'SINGLE PROTEIN', b'SINGLE PROTEIN'), (b'UNCURATED', b'UNCURATED'), (b'SUBUNIT', b'SUBUNIT')])),
                ('stoichiometry', chembl_core_db.db.customFields.ChemblPositiveIntegerField(blank=True, null=True, choices=[(0, b'0'), (1, b'1'), (2, b'2'), (3, b'3'), (12, b'12')])),
                ('targcomp_id', chembl_core_db.db.customFields.ChemblAutoField(help_text='Primary key.', serialize=False, primary_key=True, blank=True)),
                ('homologue', chembl_core_db.db.customFields.ChemblPositiveIntegerField(default=0, help_text='Indicates that the given component is a homologue of the correct component (e.g., from a different species) when set to 1. This may be the case if the sequence for the correct protein/nucleic acid cannot be found in sequence databases.', choices=[(0, b'0'), (1, b'1'), (2, b'2')])),
                ('component', models.ForeignKey(help_text='Foreign key to the component_sequences table, indicating which components belong to the target.', to='chembl_core_model.ComponentSequences')),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'target_components'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TargetDictionary',
            fields=[
                ('tid', chembl_core_db.db.customFields.ChemblAutoField(help_text='Unique ID for the target', serialize=False, primary_key=True, blank=True)),
                ('pref_name', chembl_core_db.db.customFields.ChemblCharField(help_text='Preferred target name: manually curated', max_length=200, null=True, db_index=True, blank=True)),
                ('tax_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='NCBI taxonomy id of target', null=True, db_index=True, blank=True)),
                ('organism', chembl_core_db.db.customFields.ChemblCharField(help_text='Source organism of molecuar target or tissue, or the target organism if compound activity is reported in an organism rather than a protein or tissue', max_length=150, null=True, db_index=True, blank=True)),
                ('updated_on', chembl_core_db.db.customFields.ChemblDateField(null=True, blank=True)),
                ('updated_by', chembl_core_db.db.customFields.ChemblCharField(max_length=100, null=True, blank=True)),
                ('popularity', chembl_core_db.db.customFields.ChemblPositiveIntegerField(null=True, blank=True)),
                ('insert_date', chembl_core_db.db.customFields.ChemblDateField(default=datetime.date.today, null=True)),
                ('target_parent_type', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=100, null=True, choices=[(b'MOLECULAR', b'MOLECULAR'), (b'NON-MOLECULAR', b'NON-MOLECULAR'), (b'PROTEIN', b'PROTEIN'), (b'UNDEFINED', b'UNDEFINED')])),
                ('in_starlite', chembl_core_db.db.customFields.ChemblNullableBooleanField(default=False)),
                ('species_group_flag', chembl_core_db.db.customFields.ChemblNullableBooleanField(help_text="Flag to indicate whether the target represents a group of species, rather than an individual species (e.g., 'Bacterial DHFR'). Where set to 1, indicates that any associated target components will be a representative, rather than a comprehensive set.")),
                ('downgraded', chembl_core_db.db.customFields.ChemblNullableBooleanField(default=False)),
                ('chembl', models.ForeignKey(blank=True, to='chembl_core_model.ChemblIdLookup', help_text='ChEMBL identifier for this target (for use on web interface etc)')),
                ('component_sequences', models.ManyToManyField(to='chembl_core_model.ComponentSequences', through='chembl_core_model.TargetComponents')),
                ('docs', models.ManyToManyField(to='chembl_core_model.Docs', through='chembl_core_model.Assays')),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'target_dictionary'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TargetRelations',
            fields=[
                ('relationship', chembl_core_db.db.customFields.ChemblCharField(help_text='Relationship between two targets (e.g., SUBSET OF, SUPERSET OF, OVERLAPS WITH)', max_length=20, choices=[(b'EQUIVALENT TO', b'EQUIVALENT TO'), (b'OVERLAPS WITH', b'OVERLAPS WITH'), (b'SUBSET OF', b'SUBSET OF'), (b'SUPERSET OF', b'SUPERSET OF')])),
                ('targrel_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Primary key', serialize=False, primary_key=True)),
                ('related_target', models.ForeignKey(related_name='from', db_column=b'related_tid', to='chembl_core_model.TargetDictionary', help_text='Identifier for the target that is related to the target of interest (foreign key to target_dicitionary table)')),
                ('target', models.ForeignKey(related_name='to', db_column=b'tid', to='chembl_core_model.TargetDictionary', help_text='Identifier for target of interest (foreign key to target_dictionary table)')),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'target_relations'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TargetType',
            fields=[
                ('target_type', chembl_core_db.db.customFields.ChemblCharField(help_text='Target type (as used in target dictionary)', max_length=30, serialize=False, primary_key=True)),
                ('target_desc', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of target type', max_length=250, null=True, blank=True)),
                ('parent_type', chembl_core_db.db.customFields.ChemblCharField(help_text="Higher level classification of target_type, allowing grouping of e.g., all 'PROTEIN' targets, all 'NON-MOLECULAR' targets etc.", max_length=25, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True, 'db_table': 'target_type'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='UsanStems',
            fields=[
                ('usan_stem_id', chembl_core_db.db.customFields.ChemblPositiveIntegerField(help_text='Numeric primary key.', serialize=False, primary_key=True)),
                ('stem', chembl_core_db.db.customFields.ChemblCharField(help_text='Stem defined for use in United States Adopted Names.', max_length=100)),
                ('subgroup', chembl_core_db.db.customFields.ChemblCharField(help_text='More specific subgroup of the stem defined for use in United States Adopted Names.', max_length=100)),
                ('annotation', chembl_core_db.db.customFields.ChemblCharField(help_text='Meaning of the stem (e.g., the class of compound it applies to).', max_length=2000, null=True, blank=True)),
                ('stem_class', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=100, null=True, help_text='Indicates whether stem is used as a Prefix/Infix/Suffix.', choices=[(b'Suffix', b'Suffix'), (b'Prefix', b'Prefix'), (b'Infix', b'Infix')])),
                ('major_class', chembl_core_db.db.customFields.ChemblCharField(blank=True, max_length=100, null=True, help_text='Protein family targeted by compounds of this class (e.g., GPCR/Ion channel/Protease) where known/applicable.', choices=[(b'GPCR', b'GPCR'), (b'NR', b'NR'), (b'PDE', b'PDE'), (b'ion channel', b'ion channel'), (b'kinase', b'kinase'), (b'protease', b'protease')])),
                ('who_extra', chembl_core_db.db.customFields.ChemblNullableBooleanField(default=False, help_text='Stem not represented in USAN list, but added from WHO INN stem list (where set to 1).')),
                ('downgraded', chembl_core_db.db.customFields.ChemblNullableBooleanField(default=False, help_text='Stem no longer included in USAN listing (where set to 1).')),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'usan_stems'
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Version',
            fields=[
                ('name', chembl_core_db.db.customFields.ChemblCharField(help_text='Name of release version', max_length=20, serialize=False, primary_key=True)),
                ('creation_date', chembl_core_db.db.customFields.ChemblDateField(help_text='Date database created', null=True, blank=True)),
                ('comments', chembl_core_db.db.customFields.ChemblCharField(help_text='Description of release version', max_length=2000, null=True, blank=True)),
            ],
            options={
                'abstract': False,
                'managed': True,'db_table': 'version'
            },
            bases=(models.Model,),
        ),
        migrations.AlterUniqueTogether(
            name='usanstems',
            unique_together=set([('stem', 'subgroup')]),
        ),
        migrations.AddField(
            model_name='targetdictionary',
            name='target_type',
            field=models.ForeignKey(db_column=b'target_type', blank=True, to='chembl_core_model.TargetType', help_text='Describes whether target is a protein, an organism, a tissue etc. Foreign key to TARGET_TYPE table.', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='targetcomponents',
            name='target',
            field=models.ForeignKey(db_column=b'tid', to='chembl_core_model.TargetDictionary', help_text='Foreign key to the target_dictionary, indicating the target to which the components belong.'),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='sitecomponents',
            unique_together=set([('site', 'component', 'domain')]),
        ),
        migrations.AddField(
            model_name='researchcompanies',
            name='res_stem',
            field=models.ForeignKey(blank=True, to='chembl_core_model.ResearchStem', help_text='Foreign key to research_stem table.', null=True),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='researchcompanies',
            unique_together=set([('res_stem', 'company')]),
        ),
        migrations.AlterUniqueTogether(
            name='proteinfamilyclassification',
            unique_together=set([('l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'l8')]),
        ),
        migrations.AddField(
            model_name='moleculesynonyms',
            name='res_stem',
            field=models.ForeignKey(blank=True, to='chembl_core_model.ResearchStem', help_text='Foreign key to the research_stem table. Where a synonym is a research code, this links to further information about the company associated with that code.', null=True),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='moleculesynonyms',
            unique_together=set([('molecule', 'synonyms', 'syn_type')]),
        ),
        migrations.AddField(
            model_name='moleculedictionary',
            name='assays',
            field=models.ManyToManyField(to='chembl_core_model.Assays', null=True, through='chembl_core_model.Activities', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='moleculedictionary',
            name='chembl',
            field=models.ForeignKey(blank=True, to='chembl_core_model.ChemblIdLookup', help_text='ChEMBL identifier for this compound (for use on web interface etc)', unique=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='moleculedictionary',
            name='docs',
            field=models.ManyToManyField(to='chembl_core_model.Docs', null=True, through='chembl_core_model.CompoundRecords', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='moleculedictionary',
            name='products',
            field=models.ManyToManyField(to='chembl_core_model.Products', null=True, through='chembl_core_model.Formulations', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='moleculeatcclassification',
            name='molecule',
            field=models.ForeignKey(db_column=b'molregno', to='chembl_core_model.MoleculeDictionary', help_text='Drug to which the ATC code applies (foreign key to molecule_dictionary table)'),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='mechanismrefs',
            unique_together=set([('mechanism', 'ref_type', 'ref_id')]),
        ),
        migrations.AddField(
            model_name='journalarticles',
            name='journal',
            field=models.ForeignKey(to='chembl_core_model.Journals'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='formulations',
            name='molecule',
            field=models.ForeignKey(db_column=b'molregno', blank=True, to='chembl_core_model.MoleculeDictionary', help_text='Unique identifier of the ingredient FK to MOLECULE_DICTIONARY', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='formulations',
            name='product',
            field=models.ForeignKey(help_text='Unique identifier of the product. FK to PRODUCTS', to='chembl_core_model.Products'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='formulations',
            name='record',
            field=models.ForeignKey(help_text='Foreign key to the compound_records table.', to='chembl_core_model.CompoundRecords'),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='formulations',
            unique_together=set([('record', 'product')]),
        ),
        migrations.AddField(
            model_name='drugmechanism',
            name='molecule',
            field=models.ForeignKey(db_column=b'molregno', blank=True, to='chembl_core_model.MoleculeDictionary', help_text='Molregno for the drug (foreign key to molecule_dictionary table)', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='drugmechanism',
            name='record',
            field=models.ForeignKey(help_text='Record_id for the drug (foreign key to compound_records table)', to='chembl_core_model.CompoundRecords'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='drugmechanism',
            name='site',
            field=models.ForeignKey(blank=True, to='chembl_core_model.BindingSites', help_text='Binding site for the drug within the target (where known) - foreign key to binding_sites table', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='drugmechanism',
            name='target',
            field=models.ForeignKey(db_column=b'tid', blank=True, to='chembl_core_model.TargetDictionary', help_text='Target associated with this mechanism of action (foreign key to target_dictionary table)', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='docs',
            name='journal_id',
            field=models.ForeignKey(db_column=b'journal_id', blank=True, to='chembl_core_model.Journals', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='compoundrecords',
            name='assays',
            field=models.ManyToManyField(to='chembl_core_model.Assays', null=True, through='chembl_core_model.Activities', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='compoundrecords',
            name='doc',
            field=models.ForeignKey(help_text='Foreign key to documents table', to='chembl_core_model.Docs'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='compoundrecords',
            name='molecule',
            field=models.ForeignKey(db_column=b'molregno', blank=True, to='chembl_core_model.MoleculeDictionary', help_text='Foreign key to compounds table (compound structure)', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='compoundrecords',
            name='products',
            field=models.ManyToManyField(to='chembl_core_model.Products', null=True, through='chembl_core_model.Formulations', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='compoundrecords',
            name='src',
            field=models.ForeignKey(help_text='Foreign key to source table', to='chembl_core_model.Source'),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='componentsynonyms',
            unique_together=set([('component', 'component_synonym', 'syn_type')]),
        ),
        migrations.AddField(
            model_name='componentdomains',
            name='component',
            field=models.ForeignKey(help_text='Foreign key to the component_sequences table, indicating the molecular_component that has the given domain.', to='chembl_core_model.ComponentSequences'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='componentdomains',
            name='domain',
            field=models.ForeignKey(blank=True, to='chembl_core_model.Domains', help_text='Foreign key to the domains table, indicating the domain that is contained in the associated molecular component.', null=True),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='componentdomains',
            unique_together=set([('domain', 'component', 'start_position')]),
        ),
        migrations.AddField(
            model_name='componentclass',
            name='component',
            field=models.ForeignKey(help_text='Foreign key to component_sequences table.', to='chembl_core_model.ComponentSequences'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='componentclass',
            name='protein_class',
            field=models.ForeignKey(help_text='Foreign key to the protein_classification table.', to='chembl_core_model.ProteinClassification'),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='componentclass',
            unique_together=set([('component', 'protein_class')]),
        ),
        migrations.AlterUniqueTogether(
            name='chemblidlookup',
            unique_together=set([('entity_id', 'entity_type')]),
        ),
        migrations.AlterUniqueTogether(
            name='celldictionary',
            unique_together=set([('cell_name', 'cell_source_tax_id')]),
        ),
        migrations.AddField(
            model_name='biotherapeuticcomponents',
            name='biotherapeutics',
            field=models.ForeignKey(db_column=b'molregno', to='chembl_core_model.Biotherapeutics', help_text='Foreign key to the biotherapeutics table, indicating which biotherapeutic the component is part of.'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='biotherapeuticcomponents',
            name='component',
            field=models.ForeignKey(help_text='Foreign key to the bio_component_sequences table, indicating which component is part of the biotherapeutic.', to='chembl_core_model.BioComponentSequences'),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='biotherapeuticcomponents',
            unique_together=set([('biotherapeutics', 'component')]),
        ),
        migrations.AddField(
            model_name='bindingsites',
            name='domains',
            field=models.ManyToManyField(to='chembl_core_model.Domains', null=True, through='chembl_core_model.SiteComponents', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='bindingsites',
            name='target',
            field=models.ForeignKey(db_column=b'tid', blank=True, to='chembl_core_model.TargetDictionary', help_text='Foreign key to target_dictionary. Target on which the binding site is found.', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='atcclassification',
            name='molecules',
            field=models.ManyToManyField(to='chembl_core_model.MoleculeDictionary', through='chembl_core_model.MoleculeAtcClassification'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assays',
            name='assay_type',
            field=models.ForeignKey(db_column=b'assay_type', blank=True, to='chembl_core_model.AssayType', help_text='Assay classification, e.g. B=Binding assay, A=ADME assay, F=Functional assay', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assays',
            name='cell',
            field=models.ForeignKey(blank=True, to='chembl_core_model.CellDictionary', help_text='Foreign key to cell dictionary. The cell type or cell line used in the assay', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assays',
            name='chembl',
            field=models.ForeignKey(to='chembl_core_model.ChemblIdLookup', help_text='ChEMBL identifier for this assay (for use on web interface etc)', unique=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assays',
            name='confidence_score',
            field=models.ForeignKey(db_column=b'confidence_score', blank=True, to='chembl_core_model.ConfidenceScoreLookup', help_text='Confidence score, indicating how accurately the assigned target(s) represents the actually assay target. Foreign key to CONFIDENCE_SCORE table. 0 means uncurated/unassigned, 1 = low confidence to 9 = high confidence.', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assays',
            name='curated_by',
            field=models.ForeignKey(db_column=b'curated_by', blank=True, to='chembl_core_model.CurationLookup', help_text='Indicates the level of curation of the target assignment. Foreign key to curation_lookup table.', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assays',
            name='doc',
            field=models.ForeignKey(help_text='Foreign key to documents table', to='chembl_core_model.Docs'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assays',
            name='relationship_type',
            field=models.ForeignKey(db_column=b'relationship_type', blank=True, to='chembl_core_model.RelationshipType', help_text='Flag indicating of the relationship between the reported target in the source document and the assigned target from TARGET_DICTIONARY. Foreign key to RELATIONSHIP_TYPE table.', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assays',
            name='src',
            field=models.ForeignKey(help_text='Foreign key to source table', to='chembl_core_model.Source'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assays',
            name='target',
            field=models.ForeignKey(db_column=b'tid', blank=True, to='chembl_core_model.TargetDictionary', help_text='Target identifier to which this assay has been mapped. Foreign key to target_dictionary. From ChEMBL_15 onwards, an assay will have only a single target assigned.', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assayparameters',
            name='assay',
            field=models.ForeignKey(help_text='Foreign key to assays table. The assay to which this parameter belongs', to='chembl_core_model.Assays'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assayparameters',
            name='parameter_type',
            field=models.ForeignKey(db_column=b'parameter_type', to='chembl_core_model.ParameterType', help_text='Foreign key to parameter_type table, defining the meaning of the parameter'),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='assayparameters',
            unique_together=set([('assay', 'parameter_type')]),
        ),
        migrations.AlterUniqueTogether(
            name='activitystdslookup',
            unique_together=set([('standard_type', 'standard_units')]),
        ),
        migrations.AddField(
            model_name='activities',
            name='assay',
            field=models.ForeignKey(help_text='Foreign key to the assays table (containing the assay description)', to='chembl_core_model.Assays'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='activities',
            name='data_validity_comment',
            field=models.ForeignKey(db_column=b'data_validity_comment', blank=True, to='chembl_core_model.DataValidityLookup', help_text="Comment reflecting whether the values for this activity measurement are likely to be correct - one of 'Manually validated' (checked original paper and value is correct), 'Potential author error' (value looks incorrect but is as reported in the original paper), 'Outside typical range' (value seems too high/low to be correct e.g., negative IC50 value), 'Non standard unit type' (units look incorrect for this activity type).", null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='activities',
            name='doc',
            field=models.ForeignKey(blank=True, to='chembl_core_model.Docs', help_text='Foreign key to documents table (for quick lookup of publication details - can also link to documents through compound_records or assays table)', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='activities',
            name='molecule',
            field=models.ForeignKey(db_column=b'molregno', blank=True, to='chembl_core_model.MoleculeDictionary', help_text='Foreign key to compounds table (for quick lookup of compound structure - can also link to compounds through compound_records table)', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='activities',
            name='record',
            field=models.ForeignKey(help_text='Foreign key to the compound_records table (containing information on the compound tested)', to='chembl_core_model.CompoundRecords'),
            preserve_default=True,
        ),

        migrations.RunSQL("delete from compound_mols ;alter table compound_mols drop column ctab;alter table compound_mols add column ctab mol;"),

    ]
