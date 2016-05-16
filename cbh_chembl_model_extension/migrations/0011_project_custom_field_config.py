# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0010_auto_20150330_1309'),
    ]

    operations = [
        migrations.AddField(
            model_name='project',
            name='custom_field_config',
            field=models.ForeignKey(related_name='project', default=None, blank=True,
                                    to='cbh_chembl_model_extension.CustomFieldConfig', null=True),
            preserve_default=True,
        ),
    ]
