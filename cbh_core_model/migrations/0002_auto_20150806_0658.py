# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0001_initial'),
        ('cbh_chembl_model_extension', '0025_auto_20150807_0933'),

    ]

    operations = [
        migrations.AddField(
            model_name='projecttype',
            name='level_0',
            field=models.ForeignKey(related_name='level_0', default=None,
                                    blank=True, to='cbh_core_model.CustomFieldConfig', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='projecttype',
            name='level_1',
            field=models.ForeignKey(related_name='level_1', default=None,
                                    blank=True, to='cbh_core_model.CustomFieldConfig', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='projecttype',
            name='level_2',
            field=models.ForeignKey(related_name='level_2', default=None,
                                    blank=True, to='cbh_core_model.CustomFieldConfig', null=True),
            preserve_default=True,
        ),
    ]
