# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0032_auto_20160302_1120'),
    ]

    operations = [
        migrations.RenameField(
            model_name='cbhcompoundbatch',
            old_name='batch_number',
            new_name='project_counter',
        ),
        migrations.AlterField(
            model_name='cbhcompoundbatch',
            name='related_molregno',
            field=models.ForeignKey(related_name='batches', default=None, blank=True, to='chembl_business_model.MoleculeDictionary', null=True),
            preserve_default=True,
        ),
    ]
