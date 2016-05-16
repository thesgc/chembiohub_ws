# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('chembl_core_model', '0004_auto_20150323_0942'),
    ]

    operations = [
        migrations.AlterField(
            model_name='moleculedictionary',
            name='chembl',
            field=models.ForeignKey(null=True, blank=True, to='chembl_core_model.ChemblIdLookup', help_text='ChEMBL identifier for this compound (for use on web interface etc)', unique=True),
            preserve_default=True,
        ),
    ]
