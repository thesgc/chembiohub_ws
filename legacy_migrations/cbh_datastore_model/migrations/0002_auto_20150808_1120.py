# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django_extensions.db.fields
import django_hstore.fields
import django.utils.timezone
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0007_auto_20150808_1120'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('cbh_datastore_model', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='DataFormConfig',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('created_by', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('l0', models.ForeignKey(related_name='l0', default=None, blank=True, to='cbh_core_model.CustomFieldConfig',
                                         help_text=b'The first level in the hierarchy of data points', null=True)),
                ('l1', models.ForeignKey(related_name='l1', default=None, blank=True, to='cbh_core_model.CustomFieldConfig',
                                         help_text=b'The second level in the hierarchy of data points', null=True)),
                ('l2', models.ForeignKey(related_name='l2', default=None, blank=True, to='cbh_core_model.CustomFieldConfig',
                                         help_text=b'The third level in the hierarchy of data points', null=True)),
                ('l3', models.ForeignKey(related_name='l3', default=None,
                                         blank=True, to='cbh_core_model.CustomFieldConfig', null=True)),
                ('l4', models.ForeignKey(related_name='l4', default=None, blank=True, to='cbh_core_model.CustomFieldConfig',
                                         help_text=b'The third level in the hierarchy of data points', null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='DataPointClassificationPermission',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('created_by', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('project', models.ForeignKey(to='cbh_core_model.Project')),
            ],
            options={
                'ordering': ('-modified', '-created'),
                'abstract': False,
                'get_latest_by': 'modified',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='datapoint',
            name='custom_field_config',
            field=models.ForeignKey(
                default=None, to='cbh_core_model.CustomFieldConfig', help_text=b'The schema of this datapoint'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='datapointclassification',
            name='data_form_config',
            field=models.ForeignKey(
                default=None, to='cbh_datastore_model.DataFormConfig'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='datapointclassification',
            name='l0_permission',
            field=models.ForeignKey(
                default=None, to='cbh_datastore_model.DataPointClassificationPermission'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='datapointclassification',
            name='l3',
            field=models.ForeignKey(
                related_name='l3', default=None, blank=True, to='cbh_datastore_model.DataPoint', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='datapointclassification',
            name='l4',
            field=models.ForeignKey(
                related_name='l4', default=None, blank=True, to='cbh_datastore_model.DataPoint', null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='datapoint',
            name='project_data',
            field=django_hstore.fields.SerializedDictionaryField(
                help_text=b'Data that is related to the given custom field config is stored in project_data as JSON'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='datapoint',
            name='supplementary_data',
            field=django_hstore.fields.SerializedDictionaryField(
                help_text=b'Extra data that was uploaded but was not mapped to the project data'),
            preserve_default=True,
        ),
        migrations.AlterUniqueTogether(
            name='datapointclassification',
            unique_together=set(
                [('data_form_config', 'l0', 'l1', 'l2', 'l3', 'l4')]),
        ),
        migrations.RemoveField(
            model_name='datapointclassification',
            name='project',
        ),
    ]
