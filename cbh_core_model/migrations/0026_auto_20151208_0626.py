# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):
    """Setting the username to varchar 75 without affecting django models"""
    dependencies = [
        ('cbh_core_model', '0025_auto_20151204_0626'),
    ]

    operations = [
        migrations.RunSQL("alter table auth_user alter column username type varchar(75);")
    ]
