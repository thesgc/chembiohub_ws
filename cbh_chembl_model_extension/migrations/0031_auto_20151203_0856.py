# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0030_auto_20151111_0111'),
    ]

    operations = [
        migrations.AlterField(
            model_name='cbhcompoundmultiplebatch',
            name='uploaded_file',
            field=models.ForeignKey(default=None, blank=True, to='flowjs.FlowFile', null=True),
            preserve_default=True,
        ),
    ]
