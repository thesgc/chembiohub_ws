__author__ = 'mnowotka'

import chembl_core_model.models as core

#-----------------------------------------------------------------------------------------------------------------------

class TargetType(core.TargetType):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class TargetDictionary(core.TargetDictionary):

    haystack_index = ['pref_name']
    api_exclude = []

    def save(self, force_insert=False, force_update=False, *args, **kwargs):
        needReload = False
        if not self.pk: # if we create a brand new object
            needReload = True # we need to reload to see values inserted by trigger (https://code.djangoproject.com/ticket/901)
        try:
            super(TargetDictionary, self).save(force_insert, force_update, *args, **kwargs)
        except Exception as e:
            if 'ORA-00001' in str(e):
                raise Exception('This target previously existed in ChEMBL but was removed. '
                    'It cannot be inserted again at least not currently via curation interface.')
            else:
                raise e

        if needReload:
            from_db = self.__class__.objects.get(pk=self.pk)
            for field in self.__class__._meta.fields:
                setattr(self, field.name, getattr(from_db, field.name))

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ComponentSequences(core.ComponentSequences):

    haystack_index = ['description']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ComponentSynonyms(core.ComponentSynonyms):

    haystack_index = ['component_synonym']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class TargetComponents(core.TargetComponents):

    #haystack_index = ['relationship', 'stoichiometry']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class OrganismClass(core.OrganismClass):

    #haystack_index = ['l1', 'l2', 'l3']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ProteinFamilyClassification(core.ProteinFamilyClassification):

    haystack_index = ['protein_class_desc']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ComponentClass(core.ComponentClass):

    #api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class CellDictionary(core.CellDictionary):

    #haystack_index = ['cell_name', 'cell_description', 'cell_source_tissue', 'cell_source_organism']
    api_exclude = []

    class Meta:
        proxy = True

#-----------------------------------------------------------------------------------------------------------------------

class TargetRelations(core.TargetRelations):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ProteinClassification(core.ProteinClassification):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ProteinClassSynonyms(core.ProteinClassSynonyms):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------