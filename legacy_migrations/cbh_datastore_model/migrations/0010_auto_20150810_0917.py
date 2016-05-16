# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import django_extensions.db.fields


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0010_auto_20150809_0021'),
        ('cbh_datastore_model', '0009_auto_20150810_0915'),
    ]

    operations = [
        migrations.CreateModel(
            name='DataPointClassificationPermission',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('data_point_classification', models.ForeignKey(
                    related_name='l0_permission', to='cbh_datastore_model.DataPointClassification')),
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
            model_name='datapointclassification',
            name='l0_permitted_projects',
            field=models.ManyToManyField(
                to='cbh_core_model.Project', through='cbh_datastore_model.DataPointClassificationPermission'),
            preserve_default=True,
        ),
    ]
