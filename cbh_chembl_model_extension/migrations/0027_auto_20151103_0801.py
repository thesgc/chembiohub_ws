# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0026_auto_20151103_0749'),
    ]

    operations = [
        migrations.RenameField(
            model_name='cbhcompoundbatch',
            old_name='big_image',
            new_name='bigimage',
        ),
    ]
