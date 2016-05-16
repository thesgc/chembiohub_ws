# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django_hstore.fields


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model', '0013_query'),
    ]

    operations = [
        migrations.AddField(
            model_name='query',
            name='filter',
            field=django_hstore.fields.SerializedDictionaryField(default={}),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='query',
            name='aggs',
            field=django_hstore.fields.SerializedDictionaryField(default={}),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='query',
            name='query',
            field=django_hstore.fields.SerializedDictionaryField(default={}),
            preserve_default=True,
        ),
    ]
