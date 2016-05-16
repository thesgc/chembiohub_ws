# -*- coding: utf-8 -*-
from django.db import models
from django_extensions.db.models import TimeStampedModel
from django_hstore import hstore



class Attachment(TimeStampedModel):
    number_of_rows = models.IntegerField(default=0)
    flowfile = models.ForeignKey(
        "flowjs.FlowFile", null=True, blank=True, default=None)
    created_by = models.ForeignKey("auth.User")
    data_point_classification = models.ForeignKey(
        "cbh_datastore_model.DataPointClassification")
    chosen_data_form_config = models.ForeignKey("cbh_core_model.DataFormConfig",
                                                help_text="The template data form config whose last level corresponds to the data being added", null=True, default=None, blank=True)
    attachment_custom_field_config = models.ForeignKey("cbh_core_model.CustomFieldConfig",
                                                       help_text="The schema of the table in this attachment", null=True, default=None, blank=True)
    sheet_name = models.CharField(max_length=100, default="")


class DataPoint(TimeStampedModel):

    '''Holds data for one hierachical section of a dataset. 
    Data that is related to the given custom field config is stored in project_data'''
    custom_field_config = models.ForeignKey("cbh_core_model.CustomFieldConfig",
                                            help_text="The schema of this datapoint")
    created_by = models.ForeignKey("auth.User")
    project_data = hstore.SerializedDictionaryField(
        help_text="Data that is related to the given custom field config is stored in project_data as JSON")
    supplementary_data = hstore.SerializedDictionaryField(
        help_text="Extra data that was uploaded but was not mapped to the project data")


class DataPointClassification(TimeStampedModel):
    created_by = models.ForeignKey("auth.User")
    description = models.CharField(
        max_length=1000, null=True, blank=True, default=None)
    data_form_config = models.ForeignKey("cbh_core_model.DataFormConfig")
    l0 = models.ForeignKey(DataPoint, related_name='l0', default=1)
    l1 = models.ForeignKey(DataPoint, related_name='l1', default=1)
    l2 = models.ForeignKey(DataPoint, related_name='l2', default=1)
    l3 = models.ForeignKey(DataPoint, related_name='l3', default=1)
    l4 = models.ForeignKey(DataPoint, related_name='l4', default=1)
    l0_permitted_projects = models.ManyToManyField(
        "cbh_core_model.Project", through='cbh_datastore_model.DataPointClassificationPermission')
    parent = models.ForeignKey(
        "self", null=True, default=None, related_name="children")

    class Meta:
        unique_together = (('data_form_config', 'l0', 'l1', 'l2', 'l3', 'l4'))
        ordering = ['-modified']

    def level_from(self):
        level_from = ""
        if self.l4_id != 1:
            return "l4"
        if self.l3_id != 1:
            return "l3"
        if self.l2_id != 1:
            return "l2"
        if self.l1_id != 1:
            return "l1"
        if self.l0_id != 1:
            return "l0"
        return level_from

    def all_child_generations_filter_dict(self):
        levels = ["l0", "l1", "l2", "l3", "l4"]
        levels = ["%s_id" % l for l in levels]
        used_levels = {}
        for lev in levels:
            if getattr(self, lev) != 1:
                used_levels[lev] = getattr(self, lev)
        return used_levels


class DataPointClassificationPermission(TimeStampedModel):
    project = models.ForeignKey("cbh_core_model.Project")
    data_point_classification = models.ForeignKey(
        "cbh_datastore_model.DataPointClassification", related_name="l0_permission")


class Query(TimeStampedModel):
    created_by = models.ForeignKey("auth.User")
    query = hstore.SerializedDictionaryField(default={})
    filter = hstore.SerializedDictionaryField(default={})
    aggs = hstore.SerializedDictionaryField(default={})


# def add_parent_to_datapoint_classification(sender, instance, created, **kwargs):
#     '''After saving the project make sure it has an l0 datapoint classification with th'''
#     if created is True:
#         instance.sync_permissions()
#         instance.make_editor(instance.created_by)

# post_save.connect(sync_permissions, sender=Project, dispatch_uid="proj_perms")
