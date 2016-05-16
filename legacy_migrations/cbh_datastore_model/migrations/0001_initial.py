# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django_hstore.fields
import django_extensions.db.fields
import django.utils.timezone
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0002_auto_20150806_0658'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='DataPoint',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('project_data',
                 django_hstore.fields.SerializedDictionaryField()),
                ('supplementary_data',
                 django_hstore.fields.SerializedDictionaryField()),
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
            name='DataPointClassification',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('description', models.CharField(
                    default=None, max_length=1000, null=True, blank=True)),
                ('created_by', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('l0', models.ForeignKey(related_name='l0', default=None,
                                         blank=True, to='cbh_datastore_model.DataPoint', null=True)),
                ('l1', models.ForeignKey(related_name='l1', default=None,
                                         blank=True, to='cbh_datastore_model.DataPoint', null=True)),
                ('l2', models.ForeignKey(related_name='l2', default=None,
                                         blank=True, to='cbh_datastore_model.DataPoint', null=True)),
                ('project', models.ForeignKey(to='cbh_core_model.Project')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AlterUniqueTogether(
            name='datapointclassification',
            unique_together=set([('project', 'l0', 'l1', 'l2')]),
        ),
    ]
