# -*- coding: utf-8 -*-

from django.core.management.base import BaseCommand, CommandError
try:
    # django >= 1.7
    from django.apps import apps
    get_model = apps.get_model
except ImportError:
    # django < 1.7
    from django.db.models import get_model

#-------------------------------------------------------------------------


dataset = [["AC50", "nM", "Concentration required for 50% activity"],
           ["ALB", "ug.mL-1",
            "A measurement of the albumin protein in a biological specimen."],
           ["ALBGLOB", "-",
               "The ratio of albumin to globulin in a biological specimen."],
           ["ALP", "U.L-1,IU.L-1",
            "A measurement of the alkaline phosphatase in a biological specimen."],
           ["ALT", "U.L-1,IU.L-1",
            "A measurement of the alanine aminotransferase in a biological specimen."],
           ["APTT", "s",
            "A measurement of the length of time that it takes for clotting to occur when reagents are added to a plasma specimen. The test is partial due to the absence of tissue factor (Factor III) from the reaction mixture."],
           ["AST", "U.L-1,IU.L-1",
            "A measurement of the aspartate aminotransferase in a biological specimen."],
           ["AUC", "uM.hr,ng.hr.mL-1",
            "Area under the drug concentration time curve"],
           ["BASO", "cells.uL-1",
            "A measurement of the basophils in a biological specimen."],
           ["BASOLE", "%",
            "A relative measurement (ratio or percentage) of the basophils to leukocytes in a biological specimen."],
           ["BILDIR", "ug.mL-1",
            "A measurement of the conjugated or water-soluble bilirubin in a biological specimen."],
           ["BILI", "ug.mL-1",
            "A measurement of the total bilirubin in a biological specimen."],
           ["BUN", "ug.mL-1",
            "A measurement of the urea nitrogen in a blood specimen."],
           ["CALCIUM", "ug.mL-1",
            "A measurement of the calcium in a biological specimen."],
           ["CC50", "nM,ug.mL-1",
               "Concentration required for 50% cytotoxicity"],
           ["CHLORIDE", "mEq.L-1",
            "A measurement of the chloride in a biological specimen."],
           ["CHOL", "ug.mL-1",
            "A measurement of the cholesterol in a biological specimen."],
           ["CK", "U.L-1",
            "A measurement of the total creatine kinase in a biological specimen."],
           ["CL",
               "mL.min-1.g-1,mL.min-1.kg-1,uL.min-1.(10^6cells)-1", "Clearance"],
           ["CL_renal", "mL.min-1.kg-1", "Renal clearance"],
           ["CL/F", "mL.min-1.kg-1", "Clearance/Bioavailability"],
           ["CLogD", "-", "Calculated Log(Distribution Coefficient)"],
           ["CLogD7.4", "-",
               "Calculated Log(Distribution Coefficient) at pH7.4"],
           ["CLogP", "-", "Calculated Log(Partition Coefficient)"],
           ["Cmax", "ug.mL-1,nM", "Maximum concentration"],
           ["CO2", "nM",
            "A measurement of the carbon dioxide gas in a biological specimen."],
           ["CREAT", "ug.mL-1",
            "A measurement of the creatinine in a biological specimen."],
           ["EC50", "nM,ug.mL-1", "Effective concentration for 50% activity"],
           ["ED50", "mg.kg-1,umol.kg-1", "Effective dose for 50% activity"],
           ["EOS", "cells.uL-1",
            "A measurement of the eosinophils in a biological specimen."],
           ["EOSLE", "%",
            "A relative measurement (ratio or percentage) of the eosinophils to leukocytes in a biological specimen."],
           ["F", "%", "Bioavailability"],
           ["fAUC", "uM.hr,ng.hr.mL-1",
            "Area under the free drug concentration time curve"],
           ["FIBRINO", "ug.mL-1",
            "A measurement of the fibrinogen in a biological specimen."],
           ["GGT", "IU.L-1",
            "A measurement of the gamma glutamyl transferase in a biological specimen."],
           ["GI50", "nM,ug.mL-1",
               "Concentration required for 50% growth inhibition"],
           ["GLUC", "ug.mL-1",
            "A measurement of the glucose in a biological specimen."],
           ["HCT", "%",
            "The percentage of a whole blood specimen that is composed of red blood cells (erythrocytes)."],
           ["HGB", "ug.mL-1",
            "A measurement of the hemoglobin in a biological specimen."],
           ["IC50", "nM,ug.mL-1", "Concentration required for 50% inhibition"],
           ["IC90", "nM,ug.mL-1", "Concentration required for 90% inhibition"],
           ["IC95", "nM", "Concentration required for 95% inhibition"],
           ["IC99", "nM,ug.mL-1", "Concentration required for 99% inhibition"],
           ["ID50", "mg.kg-1,umol.kg-1", "Dose required for 50% inhibition"],
           ["IFI", "%", "Inhibition Frequency (promiscuity) Index"],
           ["Inhibition", "%", "Inhibition"],
           ["IZ", "mm", "Zone of inhibition"],
           ["k_obs", "s-1", "Observed rate constant"],
           ["k_off", "s-1", "Dissociation rate constant"],
           ["k_on", "M-1.s-1", "Association rate constant"],
           ["Kd", "nM", "Dissociation constant"],
           ["Ki", "nM", "Inhibition constant"],
           ["Km", "nM", "Michaelis constant"],
           ["LC50", "nM,ug.mL-1", "Concentration required for 50% lethality"],
           ["LD50", "mg.kg-1,umol.kg-1", "Dose required for 50% lethality"],
           ["LDH", "U.L-1,IU.L-1",
            "A measurement of the lactate dehydrogenase in a biological specimen."],
           ["LIPASE", "U.L-1",
            "A measurement of the pancreatic lipase in a biological specimen."],
           ["LogP", "-", "Log(Partition Coefficient)"],
           ["LogPapp", "-", "Log(Apparent permeability)"],
           ["LYM", "cells.uL-1",
            "A measurement of the lymphocytes in a biological specimen."],
           ["LYMLE", "%",
            "A relative measurement (ratio or percentage) of the lymphocytes to leukocytes in a biological specimen."],
           ["MCC", "nM,ug.mL-1", "Minimum cytotoxic concentration"],
           ["MCH", "pg",
            "A quantitative measurement of the mean amount of hemoglobin per erythrocyte in a biological specimen."],
           ["MCHC", "%,ug.mL-1",
            "A quantitative measurement of the mean amount of hemoglobin per erythrocytes in a specified volume of a biological specimen."],
           ["MCV", "fL",
            "A quantitative measurement of the mean volume of erythrocytes in a biological specimen."],
           ["MIC", "nM,ug.mL-1", "Minimum inhibitory concentration"],
           ["MIC50", "nM,ug.mL-1",
            "Minimum inhibitory concentration in 50% of population"],
           ["MIC70", "ug.mL-1",
            "Minimum inhibitory concentration in 70% of population"],
           ["MIC80", "nM,ug.mL-1",
            "Minimum inhibitory concentration in 80% of population"],
           ["MIC90", "nM,ug.mL-1",
            "Minimum inhibitory concentration in 90% of population"],
           ["MONO", "cells.uL-1",
            "A measurement of the monocytes in a biological specimen."],
           ["MONOLE", "%",
            "A relative measurement (ratio or percentage) of the monocytes to leukocytes in a biological specimen."],
           ["MRT", "hr", "Mean residence time"],
           ["NEUTLE", "%",
            "A relative measurement (ratio or percentage) of the neutrophils to leukocytes in a biological specimen."],
           ["NEUTSG", "cells.uL-1",
            "A measurement of the segmented neutrophils in a biological specimen."],
           ["PHOS", "ug.mL-1",
            "A measurement of the phosphate in a biological specimen."],
           ["PHOSLPD", "ug.mL-1",
            "A measurement of the phospholipids in a biological specimen."],
           ["PLAT", "cells.uL-1",
            "A measurement of the platelets (non-nucleated thrombocytes) in a biological specimen."],
           ["POTASSIUM", "mEq.L-1",
            "A measurement of the potassium in a biological specimen."],
           ["Potency", "nM",
            "Concentration or dose required to elicit a specific response"],
           ["PROT", "ug.mL-1",
            "A measurement of a group of complex organic macromolecules composed of one or more alpha-L-amino acid chains in a biological specimen."],
           ["PT", "s",
            "A blood clotting measurement that evaluates the extrinsic pathway of coagulation."],
           ["RBC", "cells.uL-1",
            "A measurement of the total erythrocytes in a biological specimen."],
           ["RBCNUC", "/100WBC",
            "A measurement of the nucleated erythrocytes (large, immature nucleated erythrocytes) in a biological specimen."],
           ["RETIRBC", "%",
            "A relative measurement (ratio or percentage) of reticulocytes to erythrocytes in a biological specimen."],
           ["SODIUM", "mEq.L-1",
            "A measurement of the sodium in a biological specimen."],
           ["Solubility", "nM", "Solubility "],
           ["Solubility", "ug.mL-1", "Solubility"],
           ["T1/2", "hr", "Half-life"],
           ["TERMBW", "g",
               "The weight of a subject at a specified end point. (NCI)"],
           ["TGI", "nM,ug.mL-1",
            "Concentration required for total growth inhibition"],
           ["Tmax", "hr", "Time for maximum concentration"],
           ["TRIG", "ug.mL-1",
            "A measurement of the triglycerides in a biological specimen."],
           ["URATE", "ug.mL-1",
            "A measurement of the urate in a biological specimen."],
           ["Vd", "L.kg-1", "Volume of Distribution"],
           ["Vd/F", "L.kg-1",
            "Volume of distribution after non-intravenous administration"],
           ["Vdss", "L.kg-1", "Steady State Volume of Distribution"],
           ["Vdss/F", "L.kg-1",
            "Volume of distribution at steady state after non-intravenous administration"],
           ["WBC", "cells.uL-1",
            "A measurement of the leukocytes in a biological specimen."],
           ["WEIGHT", "g",
            "The vertical force exerted by a mass as a result of gravity. (NCI)"],
           ["XC50", "nM", "Concentration required for 50% effect"]]


from django.contrib.auth import get_user_model
from django.db import models


def insert_bioactivities():
    PinnedCustomField = get_model("cbh_core_model", "PinnedCustomField")

    CustomFieldConfig = get_model("cbh_core_model", "CustomFieldConfig")
    DataType = get_model("cbh_core_model", "DataType")
    u = get_user_model().objects.all()[0]
    dt = DataType.objects.get(name="Activity")
    for pointer in dataset:
        try:
            cfc = CustomFieldConfig.objects.create(
                name=pointer[0], created_by=u, data_type=dt)

            PinnedCustomField.objects.create(custom_field_config=cfc,
                                             name="Standard Activity Type",
                                             field_key="Standard Activity Type",
                                             position=0,
                                             description=pointer[2],
                                             default=pointer[0],
                                             field_type=PinnedCustomField.TEXT)

            PinnedCustomField.objects.create(custom_field_config=cfc,
                                             name="Standard Unit",
                                             allowed_values=pointer[1],
                                             field_key="Standard Unit",
                                             position=1,
                                             description="",
                                             default=pointer[1].split(",")[0],
                                             field_type=PinnedCustomField.UISELECT)

            PinnedCustomField.objects.create(custom_field_config=cfc,
                                             name="Standard Activity Value",
                                             allowed_values="",
                                             field_key="Standard Activity Value",
                                             position=2,
                                             description="",
                                             default="",
                                             field_type=PinnedCustomField.NUMBER)
        except:
            pass


class Command(BaseCommand):

    def handle(self, *args, **options):
        insert_bioactivities()

#-------------------------------------------------------------------------
