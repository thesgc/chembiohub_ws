# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model', '0017_attachment'),
    ]

    operations = [
        migrations.AddField(
            model_name='attachment',
            name='number_of_rows',
            field=models.IntegerField(default=0),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='attachment',
            name='chosen_data_form_config',
            field=models.ForeignKey(
                help_text=b'The template data form config whose last level corresponds to the data being added', to='cbh_core_model.DataFormConfig'),
            preserve_default=True,
        ),
    ]
