# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0029_pinnedcustomfield_visible_to_editors_only'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='pinnedcustomfield',
            name='visible_to_editors_only',
        ),
        migrations.AddField(
            model_name='pinnedcustomfield',
            name='open_or_restricted',
            field=models.CharField(default=b'open', max_length=20, choices=[(b'open', b'Open to all viewers'), (b'restricted', b'Restricted to editors')]),
            preserve_default=True,
        ),
    ]
