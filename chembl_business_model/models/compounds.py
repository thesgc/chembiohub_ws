__author__ = 'mnowotka'

from django.conf import settings
from rdkit.Chem import InchiToInchiKey
from rdkit.Chem import MolFromInchi
from rdkit.Chem import Kekulize
from rdkit.Chem import MolToMolBlock
from rdkit.Chem import MolFromMolBlock
from rdkit.Chem.rdmolfiles import MolToSmiles
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from chembl_business_model.signals import structureChanged
from chembl_business_model.exceptions import NoStandardInchi
import chembl_core_model.models as core
from raven.contrib.django.raven_compat.models import client
from chembl_business_model.utils import inchiFromPipe
from chembl_business_model.utils import cleanup
from chembl_business_model.utils import getStructure
from datetime import datetime

#-----------------------------------------------------------------------------------------------------------------------

class MoleculeDictionary(core.MoleculeDictionary):

    haystack_index = ['pref_name']
    api_exclude = []
    #defaultIndex = 'chembl_business_model'

    def save(self, force_insert=False, force_update=False, *args, **kwargs):
        needReload = False
        if not self.pk: # if we create a brand new object
            needReload = True # we need to reload to see values inserted by trigger (https://code.djangoproject.com/ticket/901)
        try:
            super(MoleculeDictionary, self).save(force_insert, force_update, *args, **kwargs)
        except Exception as e:
            if 'ORA-00001' in str(e):
                raise Exception('This compound previously existed in ChEMBL but was removed (probably it was inorganic). '
                    'This types of compounds cannot be inserted again at least not currently via curation interface.')
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

class CompoundRecords(core.CompoundRecords):

    haystack_index = ['compound_name', 'compound_key']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class RecordDrugProperties(core.RecordDrugProperties):

    haystack_index = []
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class CompoundProperties(core.CompoundProperties):

    #haystack_index = ['mw_freebase', 'alogp', 'hba', 'hbd', 'psa', 'rtb', 'ro3_pass', 'num_ro5_violations',
    #    'med_chem_friendly', 'acd_most_apka', 'acd_most_bpka', 'acd_logp', 'acd_logd', 'molecular_species', 'full_mwt']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class MoleculeHierarchy(core.MoleculeHierarchy):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ResearchStem(core.ResearchStem):

    #haystack_index = ['research_stem']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ResearchCompanies(core.ResearchCompanies):

    #haystack_index = ['company', 'country', 'previous_company']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class MoleculeSynonyms(core.MoleculeSynonyms):

    haystack_index = ['synonyms']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class Biotherapeutics(core.Biotherapeutics):

    haystack_index = ['description']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class BioComponentSequences(core.BioComponentSequences):

    haystack_index = ['description']
    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class BiotherapeuticComponents(core.BiotherapeuticComponents):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class CompoundImages(core.CompoundImages):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class CompoundMols(core.CompoundMols):

    api_exclude = []

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class CompoundStructures(core.CompoundStructures):

    #haystack_index = ['standard_inchi']
    api_exclude = []

    def save(self, force_insert=False, force_update=False, *args, **kwargs):

        changed = False
        new  =  not bool(CompoundStructures.objects.filter(pk=self.pk).count())
        if settings.OPEN_SOURCE:
            if self.molfile:
                if not new: # The structure already exists and we only want to modify it
                    super(CompoundStructures, self).save(force_insert, force_update, *args, **kwargs) # this should trigger CMPD_STR_UPDATE_TRIG, which deletes compound images and properties and nulls standard inchi, key, smiles, and molformula
                    changed = True
             #   newInchi = inchiFromPipe(self.molfile, settings.INCHI_BINARIES_LOCATION['1.02'])
                #if newInchi != self.standard_inchi:
                 #   self.standard_inchi = newInchi
                  #  changed = True
            mol = MolFromInchi(self.standard_inchi.encode("ascii"))
            if mol:
            # self.canonical_smiles = MolToSmiles(mol)
                if not self.standard_inchi:
                    raise NoStandardInchi("for CompundStructure, pk = " + str(self.pk))

                newInchiKey = InchiToInchiKey(self.standard_inchi.encode("ascii"))
                if self.standard_inchi_key != newInchiKey:
                    self.standard_inchi_key = newInchiKey
                    mol = MolFromInchi(self.standard_inchi.encode("ascii"))
                    # self.canonical_smiles = MolToSmiles(mol)
                    changed = True
                    self.molfile = MolToMolBlock(MolFromMolBlock(str(self.molfile))) # This is how we do kekulisation in RDKit...

                self.clean_fields()
                self.validate_unique()
                super(CompoundStructures, self).save(force_insert, force_update, *args, **kwargs)

        else:
            if self.molfile:
                if not new: # The structure already exists and we only want to modify it
                    super(CompoundStructures, self).save(force_insert, force_update, *args, **kwargs) # this should trigger CMPD_STR_UPDATE_TRIG, which deletes compound images and properties and nulls standard inchi, key, smiles, and molformula
                    changed = True

                data = getStructure(self.molfile)

                newInchi = data['InChI']
                if newInchi != self.standard_inchi:
                    self.standard_inchi = newInchi
                    self.standard_inchi_key = data['InChIKey']
                    #self.molformula = data['Molecular_Formula']
                    self.canonical_smiles = data['Canonical_Smiles']
                    changed = True

            if not self.standard_inchi:
                raise NoStandardInchi("for CompundStructure, pk = " + str(self.pk))

            if not self.standard_inchi_key:
                self.standard_inchi_key = InchiToInchiKey(self.standard_inchi.encode("ascii"))

            self.clean_fields()
            self.validate_unique()
            super(CompoundStructures, self).save(force_insert, force_update, *args, **kwargs)

        if changed:
            self.molecule.structure_key = self.standard_inchi_key
            self.molecule.structure_type = "MOL"
            self.molecule.molfile_update = datetime.now()
            self.molecule.save()
            structureChanged.send(sender=self.__class__, instance=self)

    class Meta:
        proxy = True
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------
