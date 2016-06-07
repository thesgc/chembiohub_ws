# -*- coding: utf-8 -*-
"""
Model classes for the chemical batches and multiple batches that are used to hold datasets in ChemBio Hub

"""

from django.db import models, connection
from django.contrib.auth.models import User
try:
    # django >= 1.7
    from django.apps import apps
    get_model = apps.get_model
except ImportError:
    # django < 1.7
    from django.db.models import get_model

#---------------
from django_hstore import hstore
import json
import random
from chembl_business_model.models import MoleculeDictionary
from chembl_business_model.models import CompoundStructures
from rdkit import Chem
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem.rdmolfiles import MolToMolBlock
from django.core.exceptions import ValidationError
from django.conf import settings
import requests
from chembl_business_model.models import CompoundStructures
from chembl_business_model.models import MoleculeDictionary
from chembl_business_model.models import ChemblIdLookup
from chembl_business_model.tasks import generateCompoundPropertiesTask
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from django_extensions.db.models import TimeStampedModel
import shortuuid
from chembl_business_model.utils import inchiFromPipe
from rdkit.Chem import InchiToInchiKey
from rdkit.Chem import MolFromInchi
from rdkit.Chem import MolToMolBlock
from picklefield.fields import PickledObjectField
from pybel import readstring
from django.db.models.signals import post_save
from cbh_chembl_model_extension.lookups import *
from copy import copy
from rdkit.Chem.AllChem import Compute2DCoords
import base64
import StringIO
from rdkit.Chem import  SDMolSupplier, MolToMolBlock, MolFromSmarts, SDMolSupplier, AllChem, Draw, SanitizeMol, SanitizeFlags,   AssignAtomChiralTagsFromStructure

from django.db import IntegrityError



hstore.DictionaryField.register_lookup(KeyValuesAny)

hstore.DictionaryField.register_lookup(KeyValuesAll)
hstore.DictionaryField.register_lookup(KeyValuesSingle)

'''
Signal code - on a webauthed system, if a user is auto-created by a webauth login,
notify a superuser.
'''
if "django_webauth" in settings.INSTALLED_APPS:
    from django.db.models.signals import post_save
    from django.core.mail import send_mail

    def email_new_user(sender, instance, **kwargs):
        """Function to email new users (oxford-specific)"""
        if kwargs["created"]:  # only for new users
            # need to email superusers to inform them a new user has logged in
            if instance.email:
                if instance.username != instance.email:
                    #This is not an invited user
                    admin_users = find_superuser()
                    email_from = 'no-reply-chemreg@chembiohub.ox.ac.uk'
                    new_user_name = '%s %s' % (
                        instance.first_name, instance.last_name)
                    for admin in admin_users:
                        html_email = '<p>Hi %s, <br>A new user has logged onto the system via Webauth, %s. <br>You should add them to any applicable projects.<br>Thanks<br>ChemBio Hub ChemReg</p>' % (
                            admin.first_name, new_user_name)
                        email_message = 'Hi %s, A new user has logged onto the system via Webauth, %s.You should add them to any applicable projects' % (
                            admin.first_name, new_user_name)
                        send_mail('New Webauth User', email_message, email_from, [
                                  admin.email, ], fail_silently=False, html_message=html_email)
                    # we also need to email the user with a welcome message
                    welcome_message = 'Welcome to ChemReg, %s! You have successfully logged in using Webauth and you will be added to projects by an admin user shortly.' % new_user_name
                    html_welcome_message = '<p>Welcome to ChemReg, %s! <br>You have successfully logged in using Webauth and you will be added to projects by an admin user shortly.</p>' % new_user_name
                    send_mail('Welcome to ChemReg', welcome_message, email_from, [
                              instance.email, ], fail_silently=False, html_message=html_welcome_message)

    def find_superuser():
        """Find the first superuser that you can possibly deprecated"""
        return User.objects.all().filter(is_active=True, is_superuser=True).exclude(email=None).exclude(email="")

    post_save.connect(email_new_user, sender=User)


def generate_uox_id():
    """Create an alphanumeric ID for a new compound or blinded compound or inventory item, 
    ensuring that this has not been created before
    There is a theoretical possibility of a race condition here but it is unlikely due to the number of possible IDs that can be generated
    We have proposed a foolproof method to do this ID generation in thecbh_chembl_id_generator repository but have not had time to implement it"""
    two_letterg = shortuuid.ShortUUID()
    two_letterg.set_alphabet("ABCDEFGHJKLMNPQRSTUVWXYZ")
    two_letter = two_letterg.random(length=2)
    two_numberg = shortuuid.ShortUUID()
    two_numberg.set_alphabet("0123456789")
    two_number = two_numberg.random(length=2)
    three_letterg = shortuuid.ShortUUID()
    three_letterg.set_alphabet("ABCDEFGHJKLMNPQRSTUVWXYZ")
    three_letter = three_letterg.random(length=3)

    uox_id = "%s%s%s%s" % (
        settings.ID_PREFIX, two_letter, two_number, three_letter)
    try:
        ChemblIdLookup.objects.get(chembl_id=uox_id)
        return generate_uox_id()
    except ObjectDoesNotExist:
        return uox_id


BONDS_WEDGED_SDF_PROP = '''>  <_drawingBondsWedged>
True

$$$$'''


#-----------------------------------------------------------------------------------------------------------------------

def _parseMolData(data):
    """Imports a molfile and verifies if all of the coordinates are set to zeros.
    if they are set to zeros then we know there are no real coordinates in the molfile
    In this case we allow RDKit to recaculate the positions of the atoms and come up with its own pictorial representation of the molecule
    If not we use the molecule as drawn"""
    suppl = SDMolSupplier()

    suppl.SetData(str(data), sanitize=False)
    data = [x for x in suppl if x]
    for x in data:
        if not x.HasProp("_drawingBondsWedged"):
            SanitizeMol(x)
        ctab = MolToMolBlock(x)
        ctablines = [item.split("0.0000") for item in ctab.split("\n") if "0.0000" in item]
        needs_redraw = 0
        for line in ctablines:
            if len(line) > 3:
                needs_redraw +=1
        if needs_redraw == len(ctablines):
             #check for overlapping molecules in the CTAB 
            SanitizeMol(x)
            Compute2DCoords(x)
    return data

def _mols2imageStream(mols, f, format, size, legend, highlightMatch=None):
    """Return an input stream for the molecule as drawn"""
    highlights = None
    if highlightMatch:
        pattern = MolFromSmarts(highlightMatch)
        highlights = [mol.GetSubstructMatch(pattern) for mol in mols]
    kek = True
    if mols[0].HasProp("_drawingBondsWedged"):
        kek=False
    fit = False
    options = DrawingOptions() 
    subim = (size,size)
    if size >150:
        subim = (size *2,size *2)
        options.coordScale = 3
        options.bondLineWidth = 3.6
        options.dblBondOffset = 0.435
        options.atomLabelFontSize = 60
        fit = True
        
    image = Draw.MolsToGridImage(mols,molsPerRow=min(len(mols),4),subImgSize=subim,
                                     kekulize=kek,highlightAtomLists=highlights, fitImage=fit,
                                    options=options
    )
    image.save(f, format)

def _mols2imageString(mols,size,legend, format, recalc=False, highlightMatch=None):
    """Take an input stream for the molecule image and return as a string"""
    if not mols:
        return ''
 #   if recalc:
  #      _apply(mols, _computeCoords)
    imageData = StringIO.StringIO()
    for mol in mols:
        try:
            SanitizeMol(mol,sanitizeOps=SanitizeFlags.SANITIZE_ALL^SanitizeFlags.SANITIZE_CLEANUPCHIRALITY^Chem.SanitizeFlags.SANITIZE_SETCONJUGATION^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        except ValueError:
            return imageData.getvalue()
        AllChem.AssignAtomChiralTagsFromStructure(mol,replaceExistingTags=False)
    _mols2imageStream(mols, imageData, format, size, legend, highlightMatch=highlightMatch)
    return imageData.getvalue()

def _ctab2image(data,size,legend, recalc=True, highlightMatch=None):
    """Generate a base64 encoded image string for the molecule"""
    data = _mols2imageString(_parseMolData(data),size,legend, 'PNG', recalc=recalc, highlightMatch=highlightMatch)
    
    return base64.b64encode(data)


def set_images(batch):
    """Append images to a particular compound batch object"""
    batch.bigimage = _ctab2image(copy(batch.ctab), 400, False, recalc=None)
    batch.image = _ctab2image(copy(batch.ctab),80,False, recalc=None)

class CBHCompoundBatchManager(models.Manager):
    """Manager methods for the generation of compound batch objects"""

    def blinded(self, project=None):
        '''Generate a batch with a blinded id'''
        return CBHCompoundBatch(project=project)

    def from_rd_mol(self, rd_mol, orig_ctab=None, smiles="", project=None, reDraw=None):
        '''Clean up the structures that come in from Smiles or from XLS or SDFs'''
        # Get a copy of the mol data
        moldata = rd_mol
        if orig_ctab is None and moldata:
            for name in moldata.GetPropNames():
                # delete the property names for the saved ctab
                moldata.ClearProp(name)
            ctab = Chem.MolToMolBlock(moldata)
            # make it into an sdf block
            ctab += "\n$$$$"
        else:
            not_first_lists = orig_ctab.split("\n")
            headerlines = True
            lines = []
            reallines = 0
            for line in not_first_lists:
                if "V2000" in line:
                    headerlines = False
                if headerlines:
                    pass
                else:
                    lines.append(line)
                    reallines += 1
                if "END" in line:
                    break
            if reallines == 2:
                raise Exception("blank_molfile")
            all_lines = ["", "", "", ] + lines + [BONDS_WEDGED_SDF_PROP, ]

            ctab = "\n".join(all_lines)

        batch = CBHCompoundBatch(ctab=ctab, original_smiles=smiles)
        set_images(batch)
        if project:
            batch.project_id = project.id
        return batch



    

def index_new_compounds():
    """Select unindexed compounds from the compound structures table and insert them into the compound mols table using the is_valid_ctab from RDKit to check if the molecule is valid for RDKit"""
    cursor = connection.cursor()
    cursor.execute(
        "INSERT INTO compound_mols (molregno , ctab) SELECT c.molregno, mol_from_ctab(molfile::cstring) ctab FROM compound_structures c LEFT OUTER JOIN compound_mols ON c.molregno = compound_mols.molregno WHERE is_valid_ctab(molfile::cstring) AND compound_mols.molregno is null;")
    return True




class CBHCompoundMultipleBatch(TimeStampedModel):
    '''Holds a list of batches'''
    created_by = models.CharField(
        max_length=50, db_index=True, null=True, blank=True, default=None, help_text="User who created this multiple batch")
    project = models.ForeignKey("cbh_core_model.Project", null=True, blank=True, default=None, help_text="Project that the multiple batch was created in")
    uploaded_data = PickledObjectField(help_text="Now used to store the headers from the file used to generate the multiple batch, picke use perhaps can be improved or is deprecated")
    uploaded_file = models.ForeignKey(
        "cbh_core_model.CBHFlowFile", null=True, blank=True, default=None, help_text="File that was uploaded to generate this multiple batch")
    batch_count = models.IntegerField(default=0)
    task_id_for_save = models.CharField(max_length=50, null=True, blank=True, default=None,)
    


class CBHCompoundBatch(TimeStampedModel):
    '''Holds the batch information for an uploaded compound before it is fully registered'''
    ctab = models.TextField(null=True, blank=True, default=None, help_text="The raw MDL molfile text (chemical table block - ctab) that was used to generate this compound batch")
    std_ctab = models.TextField(null=True, blank=True, default=None, help_text="The standardised version of the above ctab/molfile text standardised using the Inchi software")
    canonical_smiles = models.TextField(null=True, blank=True, default=None, help_text="Canonical SMILES pattern for the uploaded compound batch generated by the pybel library from openbabel")
    original_smiles = models.TextField(null=True, blank=True, default=None, help_text="The SMILES pattern that was originally uploaded")
    uncurated_fields = hstore.DictionaryField(help_text="The fields that were uploaded to the record but were not mapped on to anything, or extra information that the user wanted to include")
    image = models.TextField(default="", help_text="A small base64 encoded image for the compound, for use in the front end")
    bigimage = models.TextField(default="", help_text="A larger base64 encoded image for the compound, for use in the front end")
    created_by = models.CharField(
        max_length=50, db_index=True, null=True, blank=True, default=None, help_text="The display name for the user who created this compound batch")
    created_by_id = models.IntegerField(null=True, blank=True, default=None, help_text="The ID from the user who created this compound batch")
    related_molregno = models.ForeignKey(
        MoleculeDictionary, 
        null=True, blank=True, 
        default=None, 
        to_field="molregno", 
        help_text="The foreign key relationship to the MoleculeDictionary record generated for this compound. May be null for compounds without a structure or inventory items"
        )
    standard_inchi = models.TextField(null=True, blank=True, default=None, help_text="The standard Inchi string for the molecule as generated by the standard Inchi library")
    standard_inchi_key = models.CharField(
        max_length=50,  null=True, blank=True, default=None,  help_text="The standard Inchi key string for the molecule as generated by the standard Inchi library")
    warnings = hstore.DictionaryField(help_text="Warnings that were raised when the molecule was uploaded")
    properties = hstore.DictionaryField(help_text="System generated properties for the molecule and whether it is archived")
    custom_fields = hstore.DictionaryField(help_text="The main field where data associated with the molecule is stored - data may be strings, numbers, lists or objects")
    multiple_batch_id = models.IntegerField(default=0, help_text="The multiple batch id of this molecule and any others that were created at the same time")
    objects = CBHCompoundBatchManager()
    project = models.ForeignKey("cbh_core_model.Project", null=True, blank=True, default=None, help_text="Foreign key to the projects table, used to determine permissions on the compound batch object")
    blinded_batch_id = models.CharField(
        default="", null=True, blank=True, max_length=12, help_text="UUID assigned in cases where the compound batch is an inventory item or a molecule without a structure" )
    project_counter = models.IntegerField(default=-1, null=True, blank=True, help_text="ID for this compound batch within the context of the project")


    def generate_project_counter(self):
        """Assign an ID to the compound batch object within the given project"""
        Project = get_model("cbh_core_model", "Project")
        return Project.objects.get_next_incremental_id_for_compound(self.project_id)



    def get_uk(self,):
        """Generate a unique key for the molecule, was used in finding duplicates on upload, may now be deprecated"""
        return "%s__%d__%s" % (self.standard_inchi_key, self.project_id, "MOL")

    def save(self, *args, **kwargs):
        """Add additional generated attributes to the compound batch object as it is being saved"""
        val = kwargs.pop("validate", True)

        # self.validate()

        if self.multiple_batch_id == 0:
            mb = CBHCompoundMultipleBatch.objects.create()
            self.multiple_batch_id = mb.id
        if self.project_counter == -1:
            self.project_counter = self.generate_project_counter()
        super(CBHCompoundBatch, self).save(*args, **kwargs)
        if self.blinded_batch_id:
            uox_id_lookup = ChemblIdLookup.objects.get_or_create(chembl_id=self.blinded_batch_id,
                                                              entity_type="DOCUMENT",
                                                              entity_id=self.id)


    def validate(self, temp_props=True):
        """Could now just use standardise, deprecated"""
        self.standardise()
        set_images(self)

    def standardise(self):
        """Ensure Inchi etc was generated correctly"""
        if self.canonical_smiles:
            return
        if not self.std_ctab:
            self.std_ctab = self.ctab
        if not self.standard_inchi:

            self.standard_inchi = inchiFromPipe(
                self.std_ctab, settings.INCHI_BINARIES_LOCATION['1.02'])
        if not self.standard_inchi:
            raise Exception("inchi_error")
        else:
            self.standard_inchi_key = InchiToInchiKey(
                self.standard_inchi.encode("ascii"))

    



def generate_structure_and_dictionary(batch):
    """
    Adding the structure data to a compound batch object
    """
    print batch.__dict__
    chirality="1"
    if batch.id:
        print "not updating"
        # currently we dont update existing compound records
    else:
        if not batch.ctab:
            #blinded compound

            uox_id = generate_uox_id()
            batch.blinded_batch_id = uox_id
            batch.save(validate=False)
            
        else:
            if not batch.canonical_smiles:
                try:
                    

                    pybelmol = readstring(
                        "mol", str(batch.ctab).encode("ascii"))
                    batch.canonical_smiles = pybelmol.write(
                        "can").split("\t")[0]
                    batch.properties["cdxml"] = pybelmol.write("cdxml")
                except:
                    pass
                try:
                    mol = MolFromInchi(
                        batch.standard_inchi.encode('ascii', 'ignore'))
                    if mol:
                        batch.std_ctab = MolToMolBlock(
                            mol, includeStereo=True)
                except:
                    pass
                inchi_key = batch.standard_inchi_key
                inchi = batch.standard_inchi
                if not batch.related_molregno_id:
                    try:
                        moldict = MoleculeDictionary.objects.get(project=batch.project,
                                                                 structure_type="MOL",
                                                                 # chirality=chirality,
                                                                 structure_key=batch.standard_inchi_key)
                    except ObjectDoesNotExist:
                        print "got here"
                        uox_id = generate_uox_id()
                        rnd = random.randint(-1000000000, -2)
                        uox_id_lookup = ChemblIdLookup.objects.create(
                            chembl_id=uox_id, entity_type="COMPOUND", entity_id=rnd)

                        moldict = MoleculeDictionary.objects.get_or_create(chembl=uox_id_lookup,
                                                                           project=batch.project,
                                                                           structure_type="MOL",
                                                                           # chirality=chirality,
                                                                           structure_key=batch.standard_inchi_key)[0]
                        uox_id_lookup.entity_id = moldict.molregno
                        uox_id_lookup.save()
                        structure = CompoundStructures(
                            molecule=moldict, molfile=batch.std_ctab, standard_inchi_key=inchi_key, standard_inchi=inchi)
                        structure.save()
                        if structure.molecule_id:
                            generateCompoundPropertiesTask(structure)
                    batch.related_molregno = moldict
                batch.save(validate=False)

    return batch



