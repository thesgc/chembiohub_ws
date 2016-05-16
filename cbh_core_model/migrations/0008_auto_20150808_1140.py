# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
from django.conf import settings
import django_extensions.db.fields


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('cbh_core_model', '0007_auto_20150808_1120'),
    ]

    operations = [
        migrations.CreateModel(
            name='DataFormConfig',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('created_by', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('l0', models.ForeignKey(related_name='l0', to='cbh_core_model.CustomFieldConfig',
                                         help_text=b'The first level in the hierarchy of data points')),
                ('l1', models.ForeignKey(related_name='l1', default=None, blank=True, to='cbh_core_model.CustomFieldConfig',
                                         help_text=b'The second level in the hierarchy of data points', null=True)),
                ('l2', models.ForeignKey(related_name='l2', default=None, blank=True, to='cbh_core_model.CustomFieldConfig',
                                         help_text=b'The third level in the hierarchy of data points', null=True)),
                ('l3', models.ForeignKey(related_name='l3', default=None,
                                         blank=True, to='cbh_core_model.CustomFieldConfig', null=True)),
                ('l4', models.ForeignKey(related_name='l4', default=None, blank=True, to='cbh_core_model.CustomFieldConfig',
                                         help_text=b'The third level in the hierarchy of data points', null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AlterUniqueTogether(
            name='dataformconfig',
            unique_together=set([('l0', 'l1', 'l2', 'l3', 'l4')]),
        ),
    ]
