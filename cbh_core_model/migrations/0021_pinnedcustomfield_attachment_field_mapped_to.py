# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0020_auto_20150917_0908'),
    ]

    operations = [
        migrations.AddField(
            model_name='pinnedcustomfield',
            name='attachment_field_mapped_to',
            field=models.ForeignKey(related_name='attachment_field_mapped_from',
                                    default=None, blank=True, to='cbh_core_model.PinnedCustomField', null=True),
            preserve_default=True,
        ),
    ]
