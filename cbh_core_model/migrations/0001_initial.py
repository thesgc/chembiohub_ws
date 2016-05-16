# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
from django.conf import settings
import django_extensions.db.fields


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    state_operations = [
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
                ('schemaform', models.TextField(
                    default=b'', null=True, blank=True)),
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
                ('field_key', models.CharField(default=b'', max_length=50)),
                ('name', models.CharField(max_length=50)),
                ('description', models.CharField(
                    default=b'', max_length=1024, null=True, blank=True)),
                ('required', models.BooleanField(default=False)),
                ('part_of_blinded_key', models.BooleanField(
                    default=False, verbose_name=b'blind key')),
                ('field_type', models.CharField(default=b'char', max_length=15, choices=[(b'text', b'Short text field'), (b'number', b'Decimal field'), (b'uiselecttag', b'Choice allowing create'), (b'char', b'Short text field'), (b'imghref', b'Image link to embed'), (
                    b'href', b'Link to server or external'), (b'date', b'Date Field'), (b'integer', b'Integer field'), (b'textarea', b'Full text'), (b'uiselecttags', b'Tags field allowing create'), (b'uiselect', b'Choice field'), (b'percentage', b'Percentage field')])),
                ('allowed_values', models.CharField(
                    default=b'', max_length=1024, null=True, blank=True)),
                ('position', models.PositiveSmallIntegerField()),
                ('custom_field_config', models.ForeignKey(
                    related_name='pinned_custom_field', to='cbh_core_model.CustomFieldConfig')),
            ],
            options={
                'ordering': ['position'],
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
                ('is_default', models.BooleanField(default=False)),
                ('created_by', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('custom_field_config', models.ForeignKey(related_name='project',
                                                          default=None, blank=True, to='cbh_core_model.CustomFieldConfig', null=True)),
            ],
            options={
                'get_latest_by': 'created',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProjectType',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('name', models.CharField(
                    default=None, max_length=100, null=True, db_index=True, blank=True)),
                ('show_compounds', models.BooleanField(default=True)),
            ],
            options={
                'ordering': ('-modified', '-created'),
                'abstract': False,
                'get_latest_by': 'modified',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='SkinningConfig',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('instance_alias', models.CharField(
                    default=b'ChemiReg', max_length=50, null=True)),
                ('project_alias', models.CharField(
                    default=b'project', max_length=50, null=True)),
                ('result_alias', models.CharField(
                    default=b'result', max_length=50, null=True)),
            ],
            options={
                'verbose_name': 'Skinning Configuration',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='project',
            name='project_type',
            field=models.ForeignKey(
                default=None, blank=True, to='cbh_core_model.ProjectType', null=True),
            preserve_default=True,
        ),
    ]

    operations = [
        # By running only state operations, we are making Django think it has
        # applied this migration to the database. In reality, we renamed a
        # "cars_tires" table to "tires_tires" earlier.
        migrations.SeparateDatabaseAndState(state_operations=state_operations)
    ]
