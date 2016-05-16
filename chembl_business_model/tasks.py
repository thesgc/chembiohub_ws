__author__ = 'mnowotka'

from celery import shared_task
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Crippen
from rdkit.Chem.rdmolfiles import MolToMolBlock
from rdkit.Chem import rdMolDescriptors as Descriptors
from rdkit.Chem import Descriptors as NewDescriptors

from rdkit.Chem.SaltRemover import SaltRemover
from chembl_business_model.indigoWrapper import indigoObj
from django.core.exceptions import ValidationError
import StringIO
from django.db import IntegrityError
import requests
import base64
import decimal
from django.conf import settings
from django.db.models import Q
from chembl_business_model.utils import iterateModelRecords
from chembl_business_model.utils import iterateNModelRecords
from chembl_business_model.utils import ImageFromMolPP
from chembl_business_model.utils import cleanup
from chembl_business_model.models import CompoundImages
from chembl_business_model.models import CompoundProperties
from chembl_business_model.models import CompoundStructures
from chembl_business_model.models import CompoundRecords
from chembl_business_model.models import MoleculeHierarchy
from chembl_business_model.models import MoleculeDictionary
from chembl_business_model.models import Source
from django.db.models import Avg, Max, Min, Count
from django.db.models.fields import NOT_PROVIDED

MoleculeDictionaryDefaults = dict(map(lambda x: (x.name, x.default if x.default != NOT_PROVIDED else None), MoleculeDictionary._meta.fields))

#-----------------------------------------------------------------------------------------------------------------------

iterateModelRecordsTask = shared_task(iterateModelRecords)

#-----------------------------------------------------------------------------------------------------------------------

iterateSomeModelRecordsTask = shared_task(iterateNModelRecords)

#-----------------------------------------------------------------------------------------------------------------------

CELERY_ON = settings.CELERY_ON if hasattr(settings, 'CELERY_ON') else False

class conditional_decorator(object):
    def __init__(self, dec, condition):
        self.decorator = dec
        self.condition = condition

    def __call__(self, func):
        if not self.condition:
            # Return the function unchanged, not decorated.
            return func
        return self.decorator(func)

#-----------------------------------------------------------------------------------------------------------------------

@conditional_decorator(shared_task, CELERY_ON)
def getCompoundImageFromPipelinePilot(structure, debug=False):

    size = 500
    molfile = structure.molfile
    png_500 = base64.b64decode(ImageFromMolPP(molfile, size))

    size = 128
    png = base64.b64decode(ImageFromMolPP(molfile, size))

    molecule = structure.molecule
    if not molecule.compoundImage:
        img = CompoundImages(molecule=molecule)
    else:
        img = molecule.compoundImage

    img.png_500 = png_500
    img.png = png

    try:
        img.save()

    except IntegrityError as e:
        if debug:
            print e.message
        else:
            raise e

#-----------------------------------------------------------------------------------------------------------------------

@conditional_decorator(shared_task, CELERY_ON)
def generateCompoundImageTask(structure, debug=False):

    if debug:
        from pydev import pydevd
        pydevd.settrace('localhost', port=6901, stdoutToServer=True, stderrToServer=True)

    molecule = structure.molecule
    if not molecule.compoundImage:
        img = CompoundImages(molecule=molecule)
    else:
        img = molecule.compoundImage

    mol = Chem.MolFromMolBlock(str(structure.molfile))
    raw = Draw.MolToImage(mol, size=(500, 500))
    raw_thumb = Draw.MolToImage(mol, size=(128, 128))

    output = StringIO.StringIO()
    raw.save(output, 'PNG')
    img.png_500 = output.getvalue()
    output.close()

    output = StringIO.StringIO()
    raw_thumb.save(output, 'PNG')
    img.png = output.getvalue()
    output.close()

    try:
        img.save()

    except IntegrityError as e:
        if debug:
            print e.message
        else:
            raise e

#-----------------------------------------------------------------------------------------------------------------------

@conditional_decorator(shared_task, CELERY_ON)
def getCompoundPropertiesFromPipelinePilot(struct, debug=False):
    url = '%schembldescriptors' % settings.PIPLINE_PILOT_ENDPOINT
    result = requests.post(url, data=struct.molfile, timeout=60)
    status = result.status_code

    if status != 200:
        raise Exception("URL %s has status %s for molfile %s" % (url, status, struct.molfile))
    data = result.json()

    if not data:
        raise Exception("Empty server response for URL: %s for molfile %s" % (url, struct.molfile))

    molecule = struct.molecule
    if not molecule.compoundProperty:
        prop = CompoundProperties(molecule=molecule)
    else:
        prop = molecule.compoundProperty

    prop.mw_freebase = data.get('Molecular_Weight')
    prop.alogp = data.get('ALogP')
    prop.hba = data.get('HBA')
    prop.hbd = data.get('HBD')
    prop.psa = data.get('PSA')
    prop.rtb = data.get('RTB')
    prop.full_molformula = data.get('Molecular_Formula')
    prop.ro3_pass = data.get('RO3_PASS')
    prop.num_ro5_violations = data.get('num_ro5_violations')
    prop.med_chem_friendly = data.get('MED_CHEM_FRIENDLY')

    if data.get('ACD_MOST_APKA'):
        prop.acd_most_apka = decimal.Decimal(data['ACD_MOST_APKA'])
    else:
        prop.acd_most_apka = None

    if data.get('ACD_MOST_BPKA'):
        prop.acd_most_bpka = decimal.Decimal(data['ACD_MOST_BPKA'])
    else:
        prop.acd_most_bpka = None

    if data.get('ACD_LogP'):
        prop.acd_logp = decimal.Decimal(data['ACD_LogP'])
    else:
        prop.acd_logp = None

    if data.get('ACD_LogD'):
        prop.acd_logd = decimal.Decimal(data['ACD_LogD'])
    else:
        prop.acd_logd = None

    prop.molecular_species = data.get('Molecular_Species')
    prop.full_mwt = data.get('Molecular_Weight')

    try:
        prop.save()

    except IntegrityError as e:
        if debug:
            print e.message
        else:
            raise e

#-----------------------------------------------------------------------------------------------------------------------

@conditional_decorator(shared_task, CELERY_ON)
def generateCompoundPropertiesTask(structure, debug=False):
    if debug:
        pydevd.settrace('localhost', port=6901, stdoutToServer=True, stderrToServer=True)

    molecule = structure.molecule
    if not molecule.compoundProperty:
        prop = CompoundProperties(molecule=molecule)
    else:
        prop = molecule.compoundProperty

    saltRemover = SaltRemover()
    mol = Chem.MolFromMolBlock(str(structure.molfile))
    base = saltRemover.StripMol(mol)
    prop.hbd = Descriptors.CalcNumHBD(mol)
    prop.hba = Descriptors.CalcNumHBA(mol)
    prop.rtb = Descriptors.CalcNumRotatableBonds(mol)
    prop.alogp = Crippen.MolLogP(mol)
    prop.psa = Descriptors.CalcTPSA(mol)
    prop.full_mwt = NewDescriptors.MolWt(mol)
    # prop.exact_mass = Descriptors.CalcExactMolWt(mol)

    if base.GetNumAtoms():
        prop.mw_freebase = NewDescriptors.MolWt(base)

    
    prop.full_molformula = Descriptors.CalcMolFormula(mol)
    
    try:
        prop.save()

    except IntegrityError as e:
        if debug:
            print e.message
        else:
            raise e

#-----------------------------------------------------------------------------------------------------------------------

@conditional_decorator(shared_task, CELERY_ON)
def generateMoleculeHierarchyFromPipelinePilot(structure, debug=False):

    if debug:
        pydevd.settrace('localhost', port=6901, stdoutToServer=True, stderrToServer=True)

    molecule = structure.molecule
    if not molecule.moleculeHierarchy:
        hierarchy = MoleculeHierarchy(molecule=molecule)
    else:
        hierarchy = molecule.moleculeHierarchy

    data = cleanup(structure.molfile, 'stripsalts')

    if not data['UPDATED']:
        hierarchy.parent_molecule = molecule
    else:
        hierarchy.parent_molecule = getParentMolregnoFromBase(data['UPDATEDCTAB'])

    hierarchy.active_molecule = hierarchy.parent_molecule

    try:
        hierarchy.save()

    except IntegrityError as e:
        if debug:
            print e.message
        else:
            raise e

#-----------------------------------------------------------------------------------------------------------------------

@conditional_decorator(shared_task, CELERY_ON)
def generateMoleculeHierarchyTask(structure, debug=False):

    if debug:
        pydevd.settrace('localhost', port=6901, stdoutToServer=True, stderrToServer=True)

    molecule = structure.molecule
    if not molecule.moleculeHierarchy:
        hierarchy = MoleculeHierarchy(molecule=molecule)
    else:
        hierarchy = molecule.moleculeHierarchy

    saltRemover = SaltRemover()
    mol = Chem.MolFromMolBlock(str(structure.molfile))
    base = saltRemover.StripMol(mol)

    if mol.GetNumAtoms() == base.GetNumAtoms():
        hierarchy.parent_molecule = molecule
    else:
        hierarchy.parent_molecule = getParentMolregnoFromBase(MolToMolBlock(base))

    hierarchy.active_molecule = hierarchy.parent_molecule

    try:
        hierarchy.save()

    except IntegrityError as e:
        if debug:
            print e.message
        else:
            raise e

#-----------------------------------------------------------------------------------------------------------------------

def getParentMolregnoFromBase(molstring):

    mol = MoleculeDictionary()
    mol.save()
    struct = CompoundStructures(molecule=mol)
    struct.molfile = molstring
    try:
        struct.save()
    except ValidationError: # already exists
        mol.delete()
        if not struct.standard_inchi_key:
            raise Exception('Insane!!!')
        return CompoundStructures.objects.filter(standard_inchi_key=struct.standard_inchi_key).all()[0].molecule

    return mol

#-----------------------------------------------------------------------------------------------------------------------

def pubmedToDOI(size = 100):
    from chembl_business_model.models import Docs
    from requests import Timeout
    from requests import ConnectionError
    from django.db import transaction
    import sys

    pk = Docs._meta.pk.name
    errCount = 0
    start = 0
    timeouts = 0
    connectionErrors = 0
    count = Docs.objects.filter(doi__isnull=True).filter(pubmed_id__isnull=False).count()
    print "objects count = %s" % count

    for i in range(start, count, size):
        transaction.commit_unless_managed()
        transaction.enter_transaction_management()
        transaction.managed(True)

        if i + size > count -1:
            end = count -1
        else:
            end = i + size
        print "i =  %s, end = %s" % (i, end)
        for doc in Docs.objects.filter(doi__isnull=True).filter(pubmed_id__isnull=False).order_by(pk)[i:end]:

            try:
                res = requests.get('http://www.pmid2doi.org/rest/json/doi/%s' % doc.pubmed_id, timeout=1.0)

            except Timeout:
                timeouts += 1
                sys.stdout.write("T")
                sys.stdout.flush()
                continue

            except ConnectionError:
                doi = getDOIFromEntrez(doc.pubmed_id)
                if doi:
                    doc.doi = doi
                    doc.save()
                    sys.stdout.write("@")
                    sys.stdout.flush()
                    continue
                else:
                    connectionErrors += 1
                    sys.stdout.write("C")
                    sys.stdout.flush()
                    continue

            if res.status_code != 200:
                doi = getDOIFromEntrez(doc.pubmed_id)
                if doi:
                    doc.doi = doi
                    doc.save()
                    sys.stdout.write("@")
                    sys.stdout.flush()
                    continue
                else:
                    errCount += 1
                    print doc.pubmed_id
                    sys.stdout.flush()
                    continue

            doi = res.json()['doi']
            if len(doi) > 50:
                sys.stdout.write("L")
                sys.stdout.flush()
                continue

            doc.doi = doi
            doc.save()
            sys.stdout.write(".")
            sys.stdout.flush()

        sys.stdout.write("commit\n")
        transaction.commit()
        transaction.leave_transaction_management()
    print "finished, errors = %s, timeuouts = %s, connectionErrors = %s" % (errCount, timeouts, connectionErrors)

#-----------------------------------------------------------------------------------------------------------------------

def getDOIFromEntrez(pubmed_id):
    from Bio import Entrez
    import re

    Entrez.email = settings.ADMINS[0][1]
    handle = Entrez.efetch(db="pubmed", id=str(pubmed_id))
    text = handle.read()
    matches = re.findall(".*doi \"(\S+)\".*", text)
    if len(matches):
        return matches[0]
    return None


#-----------------------------------------------------------------------------------------------------------------------

def updateMolfileProperties(mol):

    rdp = 'recorddrugproperties__'

    aggregates = [
                    Min(rdp + 'molecule_type'),
                    #Count(rdp + 'molecule_type', distinct=True), - whta is this for???
                    Min(rdp + 'first_approval'),
                    Max(rdp + 'oral'),
                    Max(rdp + 'parenteral'),
                    Max(rdp + 'topical'),
                    Max(rdp + 'black_box_warning'),
                    Max(rdp + 'first_in_class'),
                    Max(rdp + 'chirality'),
                    Max(rdp + 'prodrug'),
                    Max(rdp + 'therapeutic_flag'),
                    Max(rdp + 'natural_product'),
                    Max(rdp + 'availability_type'),
                    Max(rdp + 'max_phase'),
                    Max(rdp + 'inorganic_flag'),
                    Max(rdp + 'usan_year'),
                    Max(rdp + 'polymer_flag'),
                ]

    sources = Source.objects.filter(src_short_name__in=('ORANGE_BOOK', 'DRUGS', 'USP/USAN'))
    salts = MoleculeHierarchy.objects\
        .filter(Q(molecule=mol) | Q(parent_molecule=mol))\
        .values_list('pk', flat=True)\
        .distinct()
    records = CompoundRecords.objects.filter(molecule__pk__in=salts).exclude(removed=1)\
        .filter(Q(src__in=sources) |
                Q(src=Source.objects.filter(src_short_name='CANDIDATES'), filename='Antibody clinical candidates'))
    res = records.aggregate(*aggregates)

    phase = res.get(rdp + 'max_phase__max',0)
    for key, value in res.items():
        attr = key.split('__')[1]
        if phase != 4:
            if hasattr(mol, attr):
                setattr(mol, attr, MoleculeDictionaryDefaults.get(attr))
        else:
            if value is None:
                continue
            if hasattr(mol, attr):
                setattr(mol, attr, value)

    mol.save()

#-----------------------------------------------------------------------------------------------------------------------

def getAltmetricScore(size = 80):
    from chembl_business_model.models import Docs
    from requests import Timeout
    from requests import ConnectionError
    from pymongo import Connection
    import sys
    import time

    connection = Connection('mongodb://mnowotka:alt5metric@ds051577.mongolab.com:51577/altmetric')
    db = connection.altmetric
    scores = db.scores
    errCount = 0
    start = 0
    notFound = 0
    timeouts = 0
    connectionErrors = 0
    count = Docs.objects.filter(doi__isnull=False).count()

    for i in range(start, count, size):
        if i + size > count -1:
            end = count -1
        else:
            end = i + size
        print "i =  %s, end = %s" % (i, end)
        for doc in Docs.objects.filter(doi__isnull=False)[i:end]:
            time.sleep(1.0)
            try:
                res = requests.get('http://api.altmetric.com/v1/doi/%s' % doc.doi, timeout=3.0, params={'key':settings.ALIMETRIC_API_KEY})
            except Timeout:
                timeouts += 1
                sys.stdout.write("T")
                sys.stdout.flush()
                continue

            except ConnectionError:
                connectionErrors += 1
                sys.stdout.write("C")
                sys.stdout.flush()
                continue

            if res.status_code != 200:
                if res.status_code == 420:
                    print "Rate limited!"
                    print "finished, errors = %s, timeuouts = %s, connectionErrors = %s, notFound = %s" % (errCount, timeouts, connectionErrors, notFound)
                    sys.exit()

                if res.status_code == 404:
                    notFound += 1
                    sys.stdout.write("N")
                    sys.stdout.flush()
                    continue

                errCount += 1
                sys.stdout.write("E")
                sys.stdout.flush()
                continue

            scores.insert(res.json())
            sys.stdout.write(".")
            sys.stdout.flush()

        sys.stdout.write("\n")
    connection.close()
    print "finished, errors = %s, timeuouts = %s, connectionErrors = %s, notFound = %s" % (errCount, timeouts, connectionErrors, notFound)

#-----------------------------------------------------------------------------------------------------------------------