# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0014_auto_20150417_1629'),
    ]

    operations = [
        migrations.AddField(
            model_name='cbhcompoundmultiplebatch',
            name='project',
            field=models.ForeignKey(
                default=None, blank=True, to='cbh_chembl_model_extension.Project', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cbhcompoundmultiplebatch',
            name='saved',
            field=models.BooleanField(default=False),
            preserve_default=True,
        ),
    ]
