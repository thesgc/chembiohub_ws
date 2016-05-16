# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import cbh_core_api.utils


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0040_flowfile_flowfilechunk'),
    ]

    operations = [
        migrations.CreateModel(
            name='CBHFlowFile',
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
                ('project', models.ForeignKey(to='cbh_core_model.Project')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CBHFlowFileChunk',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('file', models.FileField(max_length=255, upload_to=cbh_core_api.utils.chunk_upload_to)),
                ('number', models.IntegerField()),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('parent', models.ForeignKey(related_name='chunks', to='cbh_core_model.CBHFlowFile')),
            ],
            options={
                'ordering': ['number'],
            },
            bases=(models.Model,),
        ),
        migrations.RemoveField(
            model_name='flowfilechunk',
            name='parent',
        ),
        migrations.DeleteModel(
            name='FlowFile',
        ),
        migrations.DeleteModel(
            name='FlowFileChunk',
        ),
    ]
