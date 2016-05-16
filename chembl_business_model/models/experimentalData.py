__author__ = 'mnowotka'

import chembl_core_model.models as core

#-----------------------------------------------------------------------------------------------------------------------

class RelationshipType(core.RelationshipType):

    #haystack_index = ['relationship_type', 'relationship_desc']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ConfidenceScoreLookup(core.ConfidenceScoreLookup):

    #haystack_index = ['confidence_score', 'description', 'target_mapping']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class CurationLookup(core.CurationLookup):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class AssayType(core.AssayType):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class Assays(core.Assays):

    haystack_index = ['description']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class DataValidityLookup(core.DataValidityLookup):

    #haystack_index = ['activity_type', 'relation']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ParameterType(core.ParameterType):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class AssayParameters(core.AssayParameters):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class Activities(core.Activities):

#    haystack_index = ['activity_type', 'relation', 'published_value', 'published_units', 'standard_value',
#        'standard_units', 'standard_flag', 'standard_type', 'activity_comment',
#        'published_activity_type', 'manual_curation_flag', 'data_validity_comment', 'potential_duplicate']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ActivityStdsLookup(core.ActivityStdsLookup):

    #haystack_index = ['standard_type', 'definition', 'standard_units', 'normal_range_min', 'normal_range_max']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------