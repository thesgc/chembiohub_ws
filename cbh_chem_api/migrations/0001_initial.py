# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0033_auto_20160121_1126'),
    ]

    operations = [
        migrations.CreateModel(
            name='ChemregProject',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('cbh_core_model.project',),
        ),
    ]
