# -*- coding: utf-8 -*-

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
    actions = ['reindex']


    def reindex(self, request, queryset):
        from cbh_chem_api.compounds import CBHCompoundBatchResource

        cbr = CBHCompoundBatchResource()
        cbr.reindex_elasticsearch(request)
        self.message_user(request, "Successfully reindexed ChemiReg compounds")
    reindex.short_description = "Reindex all compounds"

from cbh_core_model.models import post_save, sync_permissions

post_save.connect(
    sync_permissions, sender=ChemregProject, dispatch_uid="proj_perms2")

admin.site.register(ChemregProject, ProjectAdmin)


