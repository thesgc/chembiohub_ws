# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations

def sort_compound_images(apps,stuff):
    CBHCompoundBatch = apps.get_model("cbh_chembl_model_extension", "CBHCompoundBatch")
    User = apps.get_model("auth", "User")
    from django.core.paginator import Paginator

    paginator = Paginator(CBHCompoundBatch.objects.all(), 1000) # chunks of 1000
    from cbh_chembl_model_extension.models import set_images
    for page in range(1, paginator.num_pages +1):
        for b in paginator.page(page).object_list:
            # here you can do what you want with the row
            if b.ctab:
                set_images(b)
                b.save()
        print "done page %d of %d" % (page ,paginator.num_pages)

  
class Migration(migrations.Migration):

    dependencies = [
        ('cbh_chembl_model_extension', '0029_auto_20151109_0747'),
    ]

    operations = [
        migrations.RunPython(sort_compound_images)
    ]
