# -*- coding: utf-8 -*-

from django.contrib import admin
from cbh_core_model.models import Project, PinnedCustomField, CustomFieldConfig, SkinningConfig, ProjectType, DataFormConfig

from django.contrib.admin import ModelAdmin


from django.forms.widgets import HiddenInput, TextInput
from django.db import models
from solo.admin import SingletonModelAdmin


from django import forms

class CreatedByAdmin(object):
    def save_model(self, request, obj, form, change):
        if not obj.id:
            obj.created_by = request.user
        obj.save()




class GrappelliSortableHiddenMixin(object):

    """
    Mixin which hides the sortable field with Stacked and Tabular inlines.
    This mixin must precede admin.TabularInline or admin.StackedInline.
    """
    sortable_field_name = "position"

    def formfield_for_dbfield(self, db_field, **kwargs):
        if db_field.name == self.sortable_field_name:
            kwargs["widget"] = HiddenInput()
        return super(GrappelliSortableHiddenMixin, self).formfield_for_dbfield(db_field, **kwargs)


class DataFormConfigAdmin(CreatedByAdmin, ModelAdmin):
    exclude = ["created_by", "parent"]

    def save_model(self, request, obj, form, change):
        super(DataFormConfigAdmin, self).save_model(request, obj, form, change)
        obj.get_all_ancestor_objects(request)


class PinnedCustomFieldAdmin(ModelAdmin):
    list_display = ["name",
                    "description",
                    "field_type",
                    "allowed_values", "pinned_for_datatype", "field_key"]

    exclude = ["field_key", "standardised_alias",
               "custom_field_config", "part_of_blinded_key", "position"]

    def get_queryset(self, request):
        qs = super(PinnedCustomFieldAdmin, self).get_queryset(request)
        return qs.filter(custom_field_config=None)

    def save_model(self, request, obj, form, change):
        obj.position = 0
        obj.save()


class PinnedCustomFieldInlineForm(forms.ModelForm):
    standardised_alias = forms.ModelChoiceField(required=False, queryset=PinnedCustomField.objects.exclude(
        pinned_for_datatype=None).order_by("field_key"), empty_label="Not Mapped")

    class Meta:
        model = PinnedCustomField
        exclude = [
            "field_key", "pinned_for_datatype", "attachment_field_mapped_to"]


# GrappelliSortableHiddenMixin
class PinnedCustomFieldInline(GrappelliSortableHiddenMixin, admin.TabularInline, ):
    model = PinnedCustomField
    exclude = ["field_key", "pinned_for_datatype",
               "attachment_field_mapped_to"]

    sortable_field_name = "position"
    formfield_overrides = {
        models.CharField: {'widget': TextInput(attrs={'size': '20'})},
    }
    extra = 3
    form = PinnedCustomFieldInlineForm

    def get_extra(self, request, obj=None, **kwargs):
        """Dynamically sets the number of extra forms. 0 if the related object
        already exists or the extra configuration otherwise."""
        if obj:
            # Don't add any extra forms if the related object already exists.
            return 0
        return self.extra
# Make a template have to be chosen in order to create a schema and make
# it impossible to edit schemas once created then versioning not needed


class CustomFieldConfigAdmin(ModelAdmin):

    exclude = ["created_by", ]

    search_fields = ('name',)
    ordering = ('-created',)
    date_hierarchy = 'created'
    inlines = [PinnedCustomFieldInline, ]

    def get_readonly_fields(self, request, obj=None):
        if obj:  # editing an existing object
            return self.readonly_fields + ('schemaform',)
        return self.readonly_fields

    def save_model(self, request, obj, form, change):
        obj.created_by = request.user
        obj.save()


    formfield_overrides = {
        models.CharField: {'widget': TextInput(attrs={'size': '20'})},
    }


class ProjectTypeAdmin(ModelAdmin):
    list_display = ('name', 'show_compounds')


class ProjectAdmin(CreatedByAdmin, ModelAdmin):
    prepopulated_fields = {"project_key": ("name",)}
    list_display = ('name', 'project_key', 'created', 'project_type')
    search_fields = ('name',)
    ordering = ('-created',)
    date_hierarchy = 'created'
    exclude = ["created_by"]


admin.site.register(CustomFieldConfig, CustomFieldConfigAdmin)
admin.site.register(Project, ProjectAdmin)
admin.site.register(ProjectType, ProjectTypeAdmin)
admin.site.register(SkinningConfig, SingletonModelAdmin)
admin.site.register(DataFormConfig, DataFormConfigAdmin)
admin.site.register(PinnedCustomField, PinnedCustomFieldAdmin)
