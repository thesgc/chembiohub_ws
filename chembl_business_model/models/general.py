__author__ = 'mnowotka'

import chembl_core_model.models as core

#-----------------------------------------------------------------------------------------------------------------------

class ChemblIdLookup(core.ChemblIdLookup):

    #haystack_index = ['chembl_id']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'


#-----------------------------------------------------------------------------------------------------------------------

class Version(core.Version):

    #api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------
