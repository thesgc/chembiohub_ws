# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('chembl_core_model', '0002_auto_20150323_0929'),
    ]

    operations = [
        migrations.CreateModel(
            name='DjangoCheatSheet',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('bigIntegerField', models.BigIntegerField()),
                ('booleanField', models.BooleanField()),
                ('charField', models.CharField(max_length=10)),
                ('commaSeparatedIntegerField', models.CommaSeparatedIntegerField(max_length=20)),
                ('dateField', models.DateField()),
                ('dateTimeField', models.DateTimeField()),
                ('decimalField', models.DecimalField(max_digits=9, decimal_places=3)),
                ('emailField', models.EmailField(max_length=75)),
                ('filePathField', models.FilePathField()),
                ('floatField', models.FloatField()),
                ('integerField', models.IntegerField()),
                ('ipAddressField', models.IPAddressField()),
                ('genericIPAddressField', models.GenericIPAddressField()),
                ('nullBooleanField', models.NullBooleanField()),
                ('positiveIntegerField', models.PositiveIntegerField()),
                ('positiveSmallIntegerField', models.PositiveSmallIntegerField()),
                ('slugField', models.SlugField()),
                ('smallIntegerField', models.SmallIntegerField()),
                ('textField', models.TextField()),
                ('timeField', models.TimeField()),
                ('urlField', models.URLField()),
            ],
            options={
                'abstract': False,
                'managed': True,
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ImageErrors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('error_type', models.CharField(max_length=10, choices=[(b'reg', b'Regular (big) image'), (b'tomb', b'Thumbnail')])),
                ('image', models.ForeignKey(to='chembl_core_model.CompoundImages')),
            ],
            options={
                'abstract': False,
                'managed': True,
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='InchiErrors',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('error_type', models.CharField(max_length=30, choices=[(b'1.02', b'Inchi Software version 1.02'), (b'1.03', b'Inchi Software version 1.03'), (b'1.04', b'Inchi Software version 1.04'), (b'1.04 RD', b'RDkit using Inchi ver. 1.04'), (b'1.03 RD', b'RDkit using Inchi ver. 10,3'), (b'mol read 1.4', b'Read molfile by RDkit using inchi 1.04'), (b'mol read 1.3', b'Read molfile by RDkit using inchi 1.03'), (b'open babel 2.3.2', b'Open Babel software version 2.3.2'), (b'openbabel 2.3.2 runtime', b'Open Babel 2.3.2 runtime error'), (b'indigo 1.1.5.0 linux32', b'Indigo software, version 1.1.5.0'), (b'indigo 1.1.5.0 runtime', b'Indigo software, version 1.1.5.0 runtime error')])),
                ('structure', models.ForeignKey(to='chembl_core_model.CompoundStructures')),
            ],
            options={
                'abstract': False,
                'managed': True,
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='SDF',
            fields=[
                ('originalSDF', models.FileField(upload_to=b'sdfs')),
                ('originalHash', models.CharField(unique=True, max_length=32)),
                ('cleanSDF', models.FileField(upload_to=b'sdfs')),
                ('cleanHash', models.CharField(max_length=32, serialize=False, primary_key=True)),
            ],
            options={
                'abstract': False,
                'managed': True,
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ActionType',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.actiontype',),
        ),
        migrations.CreateModel(
            name='Activities',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.activities',),
        ),
        migrations.CreateModel(
            name='ActivityStdsLookup',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.activitystdslookup',),
        ),
        migrations.CreateModel(
            name='AssayParameters',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.assayparameters',),
        ),
        migrations.CreateModel(
            name='Assays',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.assays',),
        ),
        migrations.CreateModel(
            name='AssayType',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.assaytype',),
        ),
        migrations.CreateModel(
            name='AtcClassification',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.atcclassification',),
        ),
        migrations.CreateModel(
            name='BindingSites',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.bindingsites',),
        ),
        migrations.CreateModel(
            name='BioComponentSequences',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.biocomponentsequences',),
        ),
        migrations.CreateModel(
            name='BiotherapeuticComponents',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.biotherapeuticcomponents',),
        ),
        migrations.CreateModel(
            name='Biotherapeutics',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.biotherapeutics',),
        ),
        migrations.CreateModel(
            name='CellDictionary',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.celldictionary',),
        ),
        migrations.CreateModel(
            name='ChemblIdLookup',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.chemblidlookup',),
        ),
        migrations.CreateModel(
            name='ComponentClass',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.componentclass',),
        ),
        migrations.CreateModel(
            name='ComponentDomains',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.componentdomains',),
        ),
        migrations.CreateModel(
            name='ComponentSequences',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.componentsequences',),
        ),
        migrations.CreateModel(
            name='ComponentSynonyms',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.componentsynonyms',),
        ),
        migrations.CreateModel(
            name='CompoundImages',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.compoundimages',),
        ),
        migrations.CreateModel(
            name='CompoundMols',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.compoundmols',),
        ),
        migrations.CreateModel(
            name='CompoundProperties',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.compoundproperties',),
        ),
        migrations.CreateModel(
            name='CompoundRecords',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.compoundrecords',),
        ),
        migrations.CreateModel(
            name='CompoundStructures',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.compoundstructures',),
        ),
        migrations.CreateModel(
            name='ConfidenceScoreLookup',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.confidencescorelookup',),
        ),
        migrations.CreateModel(
            name='CurationLookup',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.curationlookup',),
        ),
        migrations.CreateModel(
            name='DataValidityLookup',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.datavaliditylookup',),
        ),
        migrations.CreateModel(
            name='DefinedDailyDose',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.defineddailydose',),
        ),
        migrations.CreateModel(
            name='Docs',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.docs',),
        ),
        migrations.CreateModel(
            name='Domains',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.domains',),
        ),
        migrations.CreateModel(
            name='DrugMechanism',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.drugmechanism',),
        ),
        migrations.CreateModel(
            name='Formulations',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.formulations',),
        ),
        migrations.CreateModel(
            name='JournalArticles',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.journalarticles',),
        ),
        migrations.CreateModel(
            name='Journals',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.journals',),
        ),
        migrations.CreateModel(
            name='LigandEff',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.ligandeff',),
        ),
        migrations.CreateModel(
            name='MechanismRefs',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.mechanismrefs',),
        ),
        migrations.CreateModel(
            name='MoleculeAtcClassification',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.moleculeatcclassification',),
        ),
        migrations.CreateModel(
            name='MoleculeDictionary',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.moleculedictionary',),
        ),
        migrations.CreateModel(
            name='MoleculeHierarchy',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.moleculehierarchy',),
        ),
        migrations.CreateModel(
            name='MoleculeSynonyms',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.moleculesynonyms',),
        ),
        migrations.CreateModel(
            name='OrganismClass',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.organismclass',),
        ),
        migrations.CreateModel(
            name='ParameterType',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.parametertype',),
        ),
        migrations.CreateModel(
            name='PredictedBindingDomains',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.predictedbindingdomains',),
        ),
        migrations.CreateModel(
            name='Products',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.products',),
        ),
        migrations.CreateModel(
            name='ProteinClassification',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.proteinclassification',),
        ),
        migrations.CreateModel(
            name='ProteinClassSynonyms',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.proteinclasssynonyms',),
        ),
        migrations.CreateModel(
            name='ProteinFamilyClassification',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.proteinfamilyclassification',),
        ),
        migrations.CreateModel(
            name='RecordDrugProperties',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.recorddrugproperties',),
        ),
        migrations.CreateModel(
            name='RelationshipType',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.relationshiptype',),
        ),
        migrations.CreateModel(
            name='ResearchCompanies',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.researchcompanies',),
        ),
        migrations.CreateModel(
            name='ResearchStem',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.researchstem',),
        ),
        migrations.CreateModel(
            name='SiteComponents',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.sitecomponents',),
        ),
        migrations.CreateModel(
            name='Source',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.source',),
        ),
        migrations.CreateModel(
            name='TargetComponents',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.targetcomponents',),
        ),
        migrations.CreateModel(
            name='TargetDictionary',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.targetdictionary',),
        ),
        migrations.CreateModel(
            name='TargetRelations',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.targetrelations',),
        ),
        migrations.CreateModel(
            name='TargetType',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.targettype',),
        ),
        migrations.CreateModel(
            name='UsanStems',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.usanstems',),
        ),
        migrations.CreateModel(
            name='Version',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('chembl_core_model.version',),
        ),
    ]
