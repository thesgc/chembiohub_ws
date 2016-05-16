# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import flowjs.utils


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='FlowFile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('identifier', models.SlugField(unique=True, max_length=255)),
                ('original_filename', models.CharField(max_length=200)),
                ('total_size', models.IntegerField(default=0)),
                ('total_chunks', models.IntegerField(default=0)),
                ('total_chunks_uploaded', models.IntegerField(default=0)),
                ('state', models.IntegerField(default=1, choices=[(1, b'Uploading'), (2, b'Completed'), (3, b'Upload Error')])),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('updated', models.DateField(auto_now=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='FlowFileChunk',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('file', models.FileField(max_length=255, upload_to=flowjs.utils.chunk_upload_to)),
                ('number', models.IntegerField()),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('parent', models.ForeignKey(related_name='chunks', to='flowjs.FlowFile')),
            ],
            options={
                'ordering': ['number'],
            },
            bases=(models.Model,),
        ),
    ]
