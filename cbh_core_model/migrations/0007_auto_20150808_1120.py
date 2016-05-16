# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0006_auto_20150807_1425'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='projecttype',
            name='level_0',
        ),
        migrations.RemoveField(
            model_name='projecttype',
            name='level_1',
        ),
        migrations.RemoveField(
            model_name='projecttype',
            name='level_2',
        ),
    ]
