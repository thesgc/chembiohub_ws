from behave import given, when, then
import json

from django.db import IntegrityError
from cbh_core_ws import parser
from collections import OrderedDict


@given("I create new custom field configs and data form configs based on the data given")
def create_realdata(context):
    from cbh_core_model.models import Project, CustomFieldConfig, PinnedCustomField, ProjectType, DataType, DataFormConfig
    from cbh_datastore_model.models import DataPoint, DataPointClassification, DataPointClassificationPermission
    setup = OrderedDict(
        [("l0", {"dtype": "Project", }),
         ("l1", {"dtype": "Sub-project", }),
         ("l2", {"dtype": "Assay"}),
         ("l3", {"dtype": "Activity"}), ]
    )

    for level in setup.keys():
        setup[level]["dtypeobj"] = DataType.objects.get_or_create(
            name=setup[level]["dtype"],)
        defname = "%s data def" % setup[level]["dtype"]
        data = parser.get_custom_field_config(
             "src/cbh_datastore_ws/cbh_datastore_ws/features/fixtures/sample_data.xlsx", defname)
        setup[level]["cfc"] = CustomFieldConfig.objects.create(
            name=defname, created_by=context.user, data_type=setup[level]["dtypeobj"][0])
        for index, d in enumerate(data):
            d["custom_field_config"] = setup[level]["cfc"]
            d["position"] = index
            if d["name"] == "IC50 value":
                d["field_type"] = "decimal"
            setup[level]["cfc"].pinned_custom_field.add(
                PinnedCustomField.objects.create(**d))
        setup[level]["cfc"].save()
    df_args = {level: setup[level]["cfc"] for level in setup.keys()}
    df_args["created_by"] = context.user
    df_args["human_added"] = True
    context.dfc = DataFormConfig.objects.create(**df_args)
    context.test_case.assertEqual(context.dfc.l4_id, None)