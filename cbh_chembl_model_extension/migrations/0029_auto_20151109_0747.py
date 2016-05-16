# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations

def sort_user_ids(apps,stuff):
    CBHCompoundBatch = apps.get_model("cbh_chembl_model_extension", "CBHCompoundBatch")
    User = apps.get_model("auth", "User")
    from django.core.paginator import Paginator

    paginator = Paginator(CBHCompoundBatch.objects.all(), 1000) # chunks of 1000

    for page in range(1, paginator.num_pages +1):
        for b in paginator.page(page).object_list:
            # here you can do what you want with the row
            if b.created_by:
                if b.created_by.strip():
                    u = User.objects.filter(username=b.created_by.strip())
                    count = u.count()
                    if count >= 1:
                        b.created_by_id = u[0].id
                    else:
                        split_name = b.created_by.split(" ")
                        if len(split_name) == 2:
                            u = User.objects.filter(first_name=split_name[0].strip(), last_name=split_name[1].strip())
                            count = u.count()
                            if count == 1:
                                b.created_by_id = u[0].id
                    if b.created_by_id:
                        b.save()

  


class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0028_cbhcompoundbatch_created_by_id'),
    ]

    operations = [
        migrations.RunPython(sort_user_ids)
    ]
