# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django_hstore.fields


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0035_projecttype_saved_search_project_type'),
    ]

    operations = [
        migrations.AddField(
            model_name='projecttype',
            name='custom_field_config_template',
            field=django_hstore.fields.SerializedDictionaryField(default={}, help_text=b'A JSON template for the custom fields for a project'),
            preserve_default=True,
        ),
    ]
