# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
import django_extensions.db.fields


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0004_auto_20150807_0936'),
    ]

    operations = [
        migrations.CreateModel(
            name='DataType',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('name', models.CharField(unique=True, max_length=20)),
            ],
            options={
                'ordering': ('-modified', '-created'),
                'abstract': False,
                'get_latest_by': 'modified',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='customfieldconfig',
            name='data_type',
            field=models.ForeignKey(
                default=None, blank=True, to='cbh_core_model.DataType', null=True),
            preserve_default=True,
        ),
    ]
