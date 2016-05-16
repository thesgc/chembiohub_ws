__author__ = 'mnowotka'

import chembl_core_model.models as core

#-----------------------------------------------------------------------------------------------------------------------

class Products(core.Products):

    haystack_index = ['trade_name']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class Formulations(core.Formulations):

    #haystack_index = ['ingredient', 'strength']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class AtcClassification(core.AtcClassification):

    #haystack_index = ['who_name', 'level1', 'level2', 'level3', 'level4', 'level5', 'who_id', 'level1_description',
    #                  'level2_description', 'level3_description', 'level4_description']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class DefinedDailyDose(core.DefinedDailyDose):

    #haystack_index = ['ddd_value', 'ddd_units', 'ddd_admr', 'ddd_comment']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class UsanStems(core.UsanStems):

    #haystack_index = ['stem', 'stem_class', 'annotation', 'major_class', 'who_extra']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class MoleculeAtcClassification(core.MoleculeAtcClassification):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------