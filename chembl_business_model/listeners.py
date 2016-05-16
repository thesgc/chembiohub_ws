__author__ = 'mnowotka'

from django.dispatch import receiver
from chembl_business_model.signals import *
from chembl_business_model.models import CompoundStructures
#from chembl_business_model.models import RecordDrugProperties
from django.conf import settings
from django.db.models.signals import post_save
from tasks import generateCompoundImageTask
from tasks import generateCompoundPropertiesTask
from tasks import generateMoleculeHierarchyTask
from tasks import getCompoundImageFromPipelinePilot
from tasks import getCompoundPropertiesFromPipelinePilot
from tasks import generateMoleculeHierarchyFromPipelinePilot
from tasks import updateMolfileProperties

#-----------------------------------------------------------------------------------------------------------------------

CELERY_ON = settings.CELERY_ON if hasattr(settings, 'CELERY_ON') else False

@receiver(structureChanged, sender=CompoundStructures)
def compoundStructurePostSaveHandler(sender, **kwargs):
    struct = kwargs['instance']

    if settings.OPEN_SOURCE:
        if CELERY_ON:
            #generateCompoundImageTask.delay(struct)
            generateCompoundPropertiesTask.delay(struct)
            generateMoleculeHierarchyTask.delay(struct)

        else:
            #generateCompoundImageTask(struct)
            generateCompoundPropertiesTask(struct)
            generateMoleculeHierarchyTask(struct)
    else:
        if CELERY_ON:
            #getCompoundImageFromPipelinePilot.delay(struct)
            getCompoundPropertiesFromPipelinePilot.delay(struct)
            generateMoleculeHierarchyFromPipelinePilot.delay(struct)
        else:
            #getCompoundImageFromPipelinePilot(struct)
            getCompoundPropertiesFromPipelinePilot(struct)
            generateMoleculeHierarchyFromPipelinePilot(struct)

#-----------------------------------------------------------------------------------------------------------------------

#@receiver(post_save, sender=RecordDrugProperties)
#def recordDrugPropertiesPostSaveHandler(sender, **kwargs):
#    prop = kwargs['instance']
#    updateMolfileProperties(prop.record.molecule)

#-----------------------------------------------------------------------------------------------------------------------
