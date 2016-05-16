__author__ = 'mnowotka'

import chembl_core_model.models as core

#-----------------------------------------------------------------------------------------------------------------------

class Journals(core.Journals):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class JournalArticles(core.JournalArticles):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class Docs(core.Docs):

    haystack_index = ['title', 'abstract']
    api_exclude = []

    def save(self, force_insert=False, force_update=False, *args, **kwargs):
        needReload = False
        if not self.pk: # if we create a brand new object
            needReload = True # we need to reload to see values inserted by trigger (https://code.djangoproject.com/ticket/901)
        try:
            super(Docs, self).save(force_insert, force_update, *args, **kwargs)
        except Exception as e:
                raise e

        if needReload:
            from_db = self.__class__.objects.get(pk=self.pk)
            for field in self.__class__._meta.fields:
                setattr(self, field.name, getattr(from_db, field.name))

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class Source(core.Source):

    #haystack_index = ['src_description', 'src_short_name']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------