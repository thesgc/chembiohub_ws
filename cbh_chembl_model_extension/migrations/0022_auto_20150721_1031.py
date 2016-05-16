# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import django_extensions.db.fields


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0021_auto_20150721_0551'),
    ]

    operations = [
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
        migrations.AddField(
            model_name='project',
            name='project_type',
            field=models.ForeignKey(related_name='project', default=None,
                                    blank=True, to='cbh_chembl_model_extension.ProjectType', null=True),
            preserve_default=True,
        ),
    ]
