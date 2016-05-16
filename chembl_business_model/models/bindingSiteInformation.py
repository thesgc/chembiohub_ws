__author__ = 'mnowotka'

import chembl_core_model.models as core

#-----------------------------------------------------------------------------------------------------------------------

class BindingSites(core.BindingSites):

    haystack_index = ['site_name']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class Domains(core.Domains):

    #haystack_index = ['domain_id', 'domain_type', 'source_domain_id', 'domain_name', 'domain_description']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ComponentDomains(core.ComponentDomains):

    #haystack_index = ['start_position', 'end_position']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class SiteComponents(core.SiteComponents):

    #haystack_index = ['sitecomp_id', 'site_residues']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------