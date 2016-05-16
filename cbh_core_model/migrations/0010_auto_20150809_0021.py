# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_core_model', '0009_project_enabled_forms'),
    ]

    operations = [
        migrations.AlterField(
            model_name='dataformconfig',
            name='l0',
            field=models.ForeignKey(related_name='l0', to='cbh_core_model.CustomFieldConfig',
                                    help_text=b'The first level in the hierarchy of the form you are trying to create. For example, if curating industries, companies,  employees , teams and departments, l0 would be industries.'),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='dataformconfig',
            name='l1',
            field=models.ForeignKey(related_name='l1', default=None, blank=True, to='cbh_core_model.CustomFieldConfig',
                                    help_text=b'The second level in the hierarchy of the form you are trying to create. For example, if curating industries, companies,  employees , teams and departments, l1 would be companies.', null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='dataformconfig',
            name='l2',
            field=models.ForeignKey(related_name='l2', default=None, blank=True, to='cbh_core_model.CustomFieldConfig',
                                    help_text=b'The third level in the hierarchy of the form you are trying to create. For example, if curating industries, companies,  employees , teams and departments, l2 would be departments.', null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='dataformconfig',
            name='l3',
            field=models.ForeignKey(related_name='l3', default=None, blank=True, to='cbh_core_model.CustomFieldConfig',
                                    help_text=b'The forth level in the hierarchy of the form you are trying to create. For example, if curating industries, companies,  employees , teams and departments, l3 would be teams.', null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='dataformconfig',
            name='l4',
            field=models.ForeignKey(related_name='l4', default=None, blank=True, to='cbh_core_model.CustomFieldConfig',
                                    help_text=b'The fifth level in the hierarchy of the form you are trying to create. For example, if curating industries, companies,  employees , teams and departments, l4 would be employees.', null=True),
            preserve_default=True,
        ),
    ]
