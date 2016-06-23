# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django_extensions.db.fields
import picklefield.fields
import django_hstore.fields
import django.utils.timezone
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('chembl_business_model', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='CBHCompoundBatch',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('ctab', models.TextField(
                    default=None, null=True, blank=True)),
                ('std_ctab', models.TextField(
                    default=None, null=True, blank=True)),
                ('canonical_smiles', models.TextField(
                    default=None, null=True, blank=True)),
                ('original_smiles', models.TextField(
                    default=None, null=True, blank=True)),
                ('editable_by', django_hstore.fields.DictionaryField()),
                ('viewable_by', django_hstore.fields.DictionaryField()),
                ('created_by', models.CharField(
                    default=None, max_length=50, null=True, db_index=True, blank=True)),
                ('standard_inchi', models.TextField(
                    default=None, null=True, blank=True)),
                ('standard_inchi_key', models.CharField(
                    default=None, max_length=50, null=True, blank=True)),
                ('warnings', django_hstore.fields.DictionaryField()),
                ('properties', django_hstore.fields.DictionaryField()),
                ('custom_fields', django_hstore.fields.DictionaryField()),
                ('errors', django_hstore.fields.DictionaryField()),
                ('multiple_batch_id', models.IntegerField(default=0)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CBHCompoundMultipleBatch',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('created_by', models.CharField(
                    default=None, max_length=50, null=True, db_index=True, blank=True)),
                ('uploaded_data', picklefield.fields.PickledObjectField(
                    editable=False)),
                
            ],
            options={
                'ordering': ('-modified', '-created'),
                'abstract': False,
                'get_latest_by': 'modified',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CustomFieldConfig',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('name', models.CharField(unique=True, max_length=50)),
                ('created_by', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ('-modified', '-created'),
                'abstract': False,
                'get_latest_by': 'modified',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PinnedCustomField',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('name', models.CharField(max_length=50)),
                ('required', models.BooleanField(default=False)),
                ('part_of_blinded_key', models.BooleanField(default=False)),
                ('field_type', models.CharField(default=b'char', max_length=4, choices=[
                 (b'char', b'Short text field'), (b'text', b'Full text'), (b'pick', b'Choice field')])),
                ('allowed_values', models.CharField(
                    default=b'', max_length=1024, null=True, blank=True)),
                ('custom_field_config', models.ForeignKey(
                    to='cbh_chembl_model_extension.CustomFieldConfig')),
            ],
            options={
                'get_latest_by': 'created',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Project',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('name', models.CharField(
                    default=None, max_length=50, null=True, db_index=True, blank=True)),
                ('project_key', models.SlugField(
                    null=True, default=None, blank=True, unique=True)),
                ('created_by', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                # ('custom_field_config', models.OneToOneField(
                #     related_name='project', to='cbh_chembl_model_extension.CustomFieldConfig')),
            ],
            options={
                'get_latest_by': 'created',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='cbhcompoundbatch',
            name='project',
            field=models.ForeignKey(
                default=None, blank=True, to='cbh_chembl_model_extension.Project', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cbhcompoundbatch',
            name='related_molregno',
            field=models.ForeignKey(
                default=None, blank=True, to='chembl_business_model.MoleculeDictionary', null=True),
            preserve_default=True,
        ),
    ]
