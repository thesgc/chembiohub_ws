# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
from django.conf import settings
import django_extensions.db.fields


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('cbh_core_model',
         '0021_pinnedcustomfield_attachment_field_mapped_to'),
        ('flowjs', '0001_initial'),
        ('cbh_datastore_model', '0016_auto_20150915_0757'),
    ]

    operations = [
        migrations.CreateModel(
            name='Attachment',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(
                    default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('sheet_name', models.CharField(default=b'', max_length=100)),
                ('attachment_custom_field_config', models.ForeignKey(
                    help_text=b'The schema of the table in this attachment', to='cbh_core_model.CustomFieldConfig')),
                ('chosen_data_form_config', models.ForeignKey(
                    to='cbh_core_model.DataFormConfig')),
                ('created_by', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('data_point_classification', models.ForeignKey(
                    to='cbh_datastore_model.DataPointClassification')),
                ('flowfile', models.ForeignKey(
                    default=None, blank=True, to='flowjs.FlowFile', null=True)),
            ],
            options={
                'ordering': ('-modified', '-created'),
                'abstract': False,
                'get_latest_by': 'modified',
            },
            bases=(models.Model,),
        ),
    ]
