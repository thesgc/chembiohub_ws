# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import django_extensions.db.fields


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0018_auto_20150521_0904'),
    ]

    operations = [
        migrations.CreateModel(
            name='SkinningConfig',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('created_by', models.CharField(
                    default=None, max_length=50, null=True, db_index=True, blank=True)),
                ('instance_alias', models.CharField(
                    default=b'ChemReg', max_length=50, null=True)),
                ('project_alias', models.CharField(
                    default=b'project', max_length=50, null=True)),
                ('result_alias', models.CharField(
                    default=b'result', max_length=50, null=True)),
            ],
            options={
                'ordering': ('-modified', '-created'),
                'abstract': False,
                'get_latest_by': 'modified',
            },
            bases=(models.Model,),
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='field_type',
            field=models.CharField(default=b'char', max_length=15, choices=[(b'char', b'Short text field'), (b'uiselect', b'Choice field'), (b'textarea', b'Full text'), (b'integer', b'Integer field'), (b'date', b'Date Field'), (
                b'text', b'Short text field'), (b'percentage', b'Percentage field'), (b'uiselecttags', b'Tags field allowing create'), (b'number', b'Decimal field'), (b'uiselecttag', b'Choice allowing create')]),
            preserve_default=True,
        ),
    ]
