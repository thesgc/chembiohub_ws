# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_datastore_model', '0015_datapointclassification_parent'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='datapointclassification',
            options={'ordering': ['-modified']},
        ),
    ]
