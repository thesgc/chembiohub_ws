# -*- coding: utf-8 -*-
"""
Admin views for ChemBio Hub Platform
This file is likely deprecated as there is no loinger a need for a separate project admin for Chemireg

"""
from django.contrib import admin
from cbh_core_model.models import Project

from django.contrib.admin import ModelAdmin

from cbh_core_api.admin import CreatedByAdmin

class ChemregProject(Project):
    
    class Meta:
        proxy = True





class ProjectAdmin(CreatedByAdmin, ModelAdmin):
    prepopulated_fields = {"project_key": ("name",)}
    list_display = ('name', 'project_key', 'created', 'project_type')
    search_fields = ('name',)
    ordering = ('-created',)
    date_hierarchy = 'created'
    exclude = ["created_by"]



from cbh_core_model.models import post_save, sync_permissions

post_save.connect(
    sync_permissions, sender=ChemregProject, dispatch_uid="proj_perms2")

admin.site.register(ChemregProject, ProjectAdmin)


