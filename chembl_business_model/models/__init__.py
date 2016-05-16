__author__ = 'mnowotka'

# Order may be important

from chembl_business_model.models.general import *
from chembl_business_model.models.meta import *
from chembl_business_model.models.sourceInformation import *
from chembl_business_model.models.compounds import *
from chembl_business_model.models.targetInformation import *
from chembl_business_model.models.bindingSiteInformation import *
from chembl_business_model.models.experimentalData import *
from chembl_business_model.models.approvedDrugData import *
from chembl_business_model.models.mechanismAnnotation import *

#register signal handlers
import chembl_business_model.listeners


