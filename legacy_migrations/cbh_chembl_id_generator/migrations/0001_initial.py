# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='CBHCompoundId',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('structure_key', models.CharField(
                    unique=True, max_length=50)),
                ('assigned_id', models.CharField(unique=True, max_length=12)),
                ('original_installation_key', models.CharField(max_length=10)),
                ('original_project_key', models.CharField(max_length=20)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
