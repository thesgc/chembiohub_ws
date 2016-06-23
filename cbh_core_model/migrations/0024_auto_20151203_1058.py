# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.utils.timezone
from django.conf import settings
import django_extensions.db.fields


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('cbh_core_model', '0023_auto_20151105_1850'),
    ]

    operations = [
        migrations.CreateModel(
            name='Invitation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', django_extensions.db.fields.CreationDateTimeField(default=django.utils.timezone.now, verbose_name='created', editable=False, blank=True)),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(default=django.utils.timezone.now, verbose_name='modified', editable=False, blank=True)),
                ('email', models.CharField(unique=True, max_length=100)),
                ('first_name', models.TextField(default=b'', null=True, blank=True)),
                ('last_name', models.TextField(default=b'', null=True, blank=True)),
                ('created_by', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('projects', models.ManyToManyField(to='cbh_core_model.Project', blank=True)),
            ],
            options={
                'ordering': ('-modified', '-created'),
                'abstract': False,
                'get_latest_by': 'modified',
            },
            bases=(models.Model,),
        ),
        migrations.AlterField(
            model_name='pinnedcustomfield',
            name='field_type',
            field=models.CharField(default=b'char', max_length=15, choices=[(b'text', b'Short text field'), (b'char', b'Short text field'), (b'textarea', b'Full text'), (b'uiselect', b'Choice field'), (b'integer', b'Integer field'), (b'number', b'Decimal field'), (b'uiselecttag', b'Choice allowing create'), (b'uiselecttags', b'Tags field allowing create'), (b'percentage', b'Percentage field'), (b'date', b'Date Field'), (b'href', b'Link to server or external'), (b'imghref', b'Image link to embed'), (b'decimal', b'Decimal field'), (b'boolean', b'checkbox'), (b'related', b'TEST')]),
            preserve_default=True,
        ),
       
    ]
