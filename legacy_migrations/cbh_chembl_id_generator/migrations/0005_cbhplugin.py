# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_id_generator', '0004_auto_20150515_0950'),
    ]

    operations = [
        migrations.CreateModel(
            name='CBHPlugin',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('full_function_name', models.CharField(max_length=100)),
                ('plugin_type', models.CharField(max_length=20, choices=[
                 (b'chemreg_on_upload', b'ChemReg (applies on upload)')])),
                ('input_json_path', models.CharField(
                    help_text=b'Based on the JSON format of a molecule produced by the ChemReg API take this item as the input argument to the plugin function', max_length=200)),
                ('output_json_path', models.CharField(
                    help_text=b'Write the output of the plugin function to this position in the JSON', max_length=200)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
