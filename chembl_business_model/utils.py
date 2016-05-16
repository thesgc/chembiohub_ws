__author__ = 'mnowotka'

from django.conf import settings
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import RDLogger
from pybel import ob, readstring
from indigoWrapper import *
import requests
from django.utils.http import urlquote
from rdkit.Chem import InchiToInchiKey
from base64 import b64encode
import hashlib

from PIL import Image
import StringIO
import os

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
ob.obErrorLog.SetOutputLevel(0)

INCHI_SPECIAL_CHARS = '={}()-/,;+?.'

#-----------------------------------------------------------------------------------------------------------------------

def getHScore():
    from django.db.models import Count
    from chembl_business_model.models import MoleculeDictionary
    from chembl_core_model.models import Assays

    ChIndex = None
    AssIndex = None

    order = MoleculeDictionary.objects.filter(downgraded=False).annotate(activities_count=Count("activities")).distinct().order_by('-activities_count')
    for idx, mol in enumerate(order):
        if mol.activities_count < idx:
            ChIndex = idx - 1
            break

    order = Assays.objects.annotate(activities_count=Count("activities")).distinct().order_by('-activities_count')
    for idx, ass in enumerate(order):
        if ass.activities_count < idx:
            AssIndex = idx - 1
            break

    return (ChIndex, AssIndex)

#-----------------------------------------------------------------------------------------------------------------------

def check_indigo_correct(size=1000):

    from chembl_business_model.models import CompoundStructures
    from clint.textui import progress
    import tempfile

    f = tempfile.NamedTemporaryFile(delete=False)
    print "saving to file %s" % f.name
    errorCount = 0
    structures = CompoundStructures.objects.all()
    count = structures.count()
    pk = CompoundStructures._meta.pk.name

    for i in progress.bar(range(0, count, size), label="Indigo check "):
        if i < 0:
            chunk = CompoundStructures.objects.order_by(pk)[:size]
        else:
            last_pk = CompoundStructures.objects.order_by(pk).only(pk).values_list(pk)[i][0]
            chunk = CompoundStructures.objects.order_by(pk).filter(pk__gt=last_pk)[:size]

        for structure in chunk:
            try:
                indigoObj.loadMolecule(str(structure.molfile))
            except Exception as e:
                f.write('%s\t%s\n' % (structure.pk, str(e)))
                errorCount += 1
    f.close()
    print "%s errors saved to %s" % (str(errorCount), f.name)

#-----------------------------------------------------------------------------------------------------------------------

def check_Activities():

    from chembl_business_model.models import Activities
    import cx_Oracle

    run = True
    excludes = [0]

    while run:
        try:
            for obj in Activities.objects.order_by(Activities._meta.pk.name).exclude(pk__lte=max(excludes)).iterator():
                pass
        except cx_Oracle.DatabaseError:
            problem = obj.pk+1
            if problem in excludes:
                problem += 1
            excludes.append(problem)
            print str(excludes)
            continue
        run = False

    print excludes

#-----------------------------------------------------------------------------------------------------------------------

def smileFromImage(image, path):
    return fromImage(image, path, smileToCanonicalSmile)

#-----------------------------------------------------------------------------------------------------------------------

def molsFromImage(image, path):
    return fromImage(image, path, smileToMol)

#-----------------------------------------------------------------------------------------------------------------------

def fromImage(image, path, fun):
    from subprocess import PIPE, Popen
    import tempfile

    fd, fpath = tempfile.mkstemp()
    os.write(fd, image)
    os.close(fd)
    arguments = [path, '-ij', '-f', 'smi', fpath]
    p = Popen(arguments, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    a, err = p.communicate(input=image)
    os.remove(fpath)

    return map(lambda x : fun(x) if x else None, filter(bool,a.split('\n')))

#-----------------------------------------------------------------------------------------------------------------------

def smileToMol(smile):
    mol = Chem.MolFromSmiles(smile)
    if mol:
        AllChem.Compute2DCoords(mol)
        return Chem.MolToMolBlock(mol)
    else:
        return None

#-----------------------------------------------------------------------------------------------------------------------

def nameToMol(name):
    res = requests.get(settings.OPSIN_URL + name + '.smi', timeout=60)
    if res.status_code != 200:
        return None
    return smileToMol(str(res.text))

#-----------------------------------------------------------------------------------------------------------------------

def smilesFromMol(mol):
    molecule = Chem.MolFromMolBlock(mol)
    smiles = Chem.MolToSmiles(molecule)
    return smiles

#-----------------------------------------------------------------------------------------------------------------------

def smileToCanonicalSmile(smile):
    mol = Chem.MolFromSmiles(smile)
    if mol:
        return Chem.MolToSmiles(mol, True)
    return smile

#-----------------------------------------------------------------------------------------------------------------------

def jsonFromSmiles(smile, size):
    mol = Chem.MolFromSmiles(smile)
    return Draw.MolToJSON(mol, size)

#-----------------------------------------------------------------------------------------------------------------------

def jsonFromMol(mol, size):
    molecule = Chem.MolFromMolBlock(mol)
    return Draw.MolToJSON(molecule, size)

#-----------------------------------------------------------------------------------------------------------------------

def inchiFromPipe(molfile, path):
    from subprocess import PIPE, Popen

    p = Popen([path, "-STDIO", "-AuxNone"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    a = p.communicate(input=str(molfile))
    return a[0][13:-1]

#-----------------------------------------------------------------------------------------------------------------------

def iterateModelRecords(modelClass, function, divFactor=300):
    all = modelClass.objects.all().count()

    for i in range(divFactor):
        low = i * (all / divFactor)
        high = (i + 1) * (all / divFactor)

        if i == (divFactor - 1):
            high = (all - 1)

        for record in modelClass.objects.all()[low:high]:
            function(record)

#-----------------------------------------------------------------------------------------------------------------------

def iterateNModelRecords(modelClass, function, N):
    for record in modelClass.objects.all()[0:N]:
        function(record)

#-----------------------------------------------------------------------------------------------------------------------

def checkInchiBinary(struct, version):
    from chembl_business_model.models import InchiErrors

    if not struct.standard_inchi or not struct.standard_inchi_key or not struct.molfile:
        return

    if struct.standard_inchi != inchiFromPipe(struct.molfile, settings.INCHI_BINARIES_LOCATION[version]):
        error = InchiErrors(error_type=version, structure=struct)
        error.save()

#-----------------------------------------------------------------------------------------------------------------------

def checkOSRA(molecule):
    img = molecule.compoundimages.png_500
    im  = Image.open(StringIO.StringIO(molecule.compoundimages.png_500))
    canonical_smiles = molecule.compoundstructures.canonical_smiles
    smile = smileFromImage(img, settings.OSRA_BINARIES_LOCATION['2.0.0'], canonical_smiles)
    im.show()
    return canonical_smiles, Chem.MolToSmiles(Chem.MolFromSmiles(smile[0]), True)

#-----------------------------------------------------------------------------------------------------------------------

def checkImage(compoundImage):
    from chembl_business_model.models import ImageErrors
    if not compoundImage.png or not compoundImage.png_500:
        return

    im = Image.open(StringIO.StringIO(compoundImage.png))
    try:
        im.verify()
        if im.size != (128, 128):
            error = ImageErrors(error_type='tomb size', image=compoundImage)
            error.save()

    except Exception as e:
        error = ImageErrors(error_type='tomb', image=compoundImage)
        error.save()

    im = Image.open(StringIO.StringIO(compoundImage.png_500))
    try:
        im.verify()
        if im.size != (500, 500):
            error = ImageErrors(error_type='reg size', image=compoundImage)
            error.save()

    except Exception as e:
        error = ImageErrors(error_type='reg', image=compoundImage)
        error.save()

#-----------------------------------------------------------------------------------------------------------------------

def getSynonymTypes():
    from chembl_business_model.models import MoleculeSynonyms
    return MoleculeSynonyms.objects.all().values_list('syn_type', flat=True).order_by('syn_type').distinct()

#-----------------------------------------------------------------------------------------------------------------------

def getImage(mol=False, molregno=False):
    from chembl_business_model.models import MoleculeDictionary
    filters = dict()
    if mol:
        key = InchiKeyFromMol(mol)
        filters['structure_key'] = key
    else:
        filters['molregno'] = molregno
    return b64encode(MoleculeDictionary.objects.filter(**filters).values_list('compoundimages__png_500')[0][0])

#-----------------------------------------------------------------------------------------------------------------------

def checkPybelInchi(struct):
    from chembl_business_model.models import InchiErrors
    if not struct.standard_inchi or not struct.standard_inchi_key or not struct.molfile:
        return

    try:
        mol = readstring('mol', str(struct.molfile))
        inchi = mol.write('inchi')

        if inchi.strip() != struct.standard_inchi:
            error = InchiErrors(error_type='open babel 2.3.2', structure=struct)
            error.save()

    except Exception:
        error = InchiErrors(error_type='openbabel 2.3.2 runtime', structure=struct)
        error.save()

#-----------------------------------------------------------------------------------------------------------------------

def checkIndigoInchi(struct):
    from chembl_business_model.models import InchiErrors
    if not struct.standard_inchi or not struct.standard_inchi_key or not struct.molfile:
        return

    try:
        mol = indigoObj.loadMolecule(str(struct.molfile))
        inchi = indigo_inchiObj.getInchi(mol)

        if inchi != struct.standard_inchi:
            error = InchiErrors(error_type='indigo 1.1.5.0 linux32', structure=struct)
            error.save()

    except Exception:
        error = InchiErrors(error_type='indigo 1.1.5.0 runtime', structure=struct)
        error.save()

#-----------------------------------------------------------------------------------------------------------------------

def checkRDkitInchi(struct):
    from chembl_business_model.models import InchiErrors
    if not struct.standard_inchi or not struct.standard_inchi_key or not struct.molfile:
        return

    m = Chem.MolFromMolBlock(str(struct.molfile))
    if not m:
        error = InchiErrors(error_type='mol read 1.4', structure=struct)
        error.save()
        return

    inchi = Chem.inchi.MolToInchi(m)
    if struct.standard_inchi != inchi:
        error = InchiErrors(error_type='1.04 RD', structure=struct)
        error.save()

#-----------------------------------------------------------------------------------------------------------------------

def tagsFromText(text):
    from sklearn.feature_extraction.text import CountVectorizer
    import numpy as np

    cv = CountVectorizer(min_df=1, charset_error="ignore",
        stop_words="english", max_features=200)
    counts = cv.fit_transform([text]).toarray().ravel()
    words = np.array(cv.get_feature_names())
    words = words[counts > 1]
    counts = counts[counts > 1]
    words = words[np.array(map(lambda x: x.isalpha(), words))]
    counts = counts[np.array(map(lambda x: x.isalpha(), words))]
    #TODO: stemming, words len > 2, remove verbs
    return [words, counts]


#-----------------------------------------------------------------------------------------------------------------------

def entitiesFromText(text):
    import jpype

    ret = set()

    jpype.startJVM(settings.JAVA_VIRTUAL_MACHINE_LOCATION, settings.OSCAR_BINARIES_LOCATION)

    Oscar = jpype.JClass("uk.ac.cam.ch.wwmm.oscar.Oscar")
    FormatType = jpype.JClass("uk.ac.cam.ch.wwmm.oscar.chemnamedict.entities.FormatType")
    oscar = Oscar()
    named_entities = oscar.findAndResolveNamedEntities(text)

    for ne in named_entities:
        smiles  = ne.getFirstChemicalStructure(FormatType.SMILES)
        if not smiles:
            continue
        smiles = smiles.getValue()
        name    = ne.getSurface()
        ne_type = ne.getType().toString()

        ret.add((name, ne_type, smiles))

    jpype.utilusM()
    return list(ret)

#-----------------------------------------------------------------------------------------------------------------------

def entitiesFromTextNew(text):

    result = requests.post(settings.OSCAR_ENDPOINT, data={'text':text, 'filter' : 'true'}, timeout=60)
    status = result.status_code

    if status != 200:
        raise Exception("URL %s has status %s" % (settings.OSCAR_ENDPOINT, status))
    return result.json()

#-----------------------------------------------------------------------------------------------------------------------

def journalChoices():
    from chembl_core_model.models import Docs
    return Docs.objects.values_list('journal', flat=True).order_by('journal').distinct()

#-----------------------------------------------------------------------------------------------------------------------

def docTypeChoices():
    from chembl_core_model.models import Docs
    return Docs.objects.values_list('doc_type', flat=True).order_by('doc_type').distinct()

#-----------------------------------------------------------------------------------------------------------------------

def metaFromDoi(doi):
    from Bio import Entrez
    from BeautifulSoup import BeautifulSoup
    from chembl_business_model.models import JournalArticles, Docs

    doc_id = None

    meta = {'journal':{'pubDate':{}}, 'authors':[]}

    Entrez.email = settings.ADMINS[0][1]
    handle = Entrez.esearch(db="pubmed", term=str(doi))
    record = BeautifulSoup(handle.read())
    id = str(record.id.getText())
    handle = Entrez.efetch(db="pubmed", id=id, rettype="gb")
    result = BeautifulSoup(handle.read())
    meta['journal']['volume'] = result.volume.getText() if result.volume else ''
    meta['journal']['issue'] = result.issue.getText() if result.issue else ''
    meta['pubmed'] = id
    meta['doi'] = result.elocationid.getText() if result.elocationid else ''
    meta['title'] = result.articletitle.getText() if result.articletitle else ''
    meta['abstract'] = result.abstracttext.getText() if result.abstracttext else ''

    journal = result.journal
    if journal:
        meta['journal']['issn'] = journal.issn.getText() if journal.issn else ''
        meta['journal']['title'] = journal.title.getText() if journal.title else ''
        meta['journal']['ISOAbbreviation'] = journal.isoabbreviation.getText() if journal.isoabbreviation else ''
        pubdate = journal.pubdate
        if pubdate:
            meta['journal']['pubDate']['year'] = pubdate.year.getText() if pubdate.year else ''
            meta['journal']['pubDate']['month'] = pubdate.month.getText() if pubdate.month else ''
            meta['journal']['pubDate']['day'] = pubdate.day.getText() if pubdate.day else ''

    if result.authorlist:
        for i in result.authorlist.childGenerator():
            if i and str(i).strip():
                author = BeautifulSoup(str(i))
                auth = {}
                if author.forename:
                    auth['forename'] = author.forename.getText()
                    auth['lastname'] = author.lastname.getText()
                    auth['initials'] = author.initials.getText()
                    meta['authors'].append(auth)

    try:
        pubmedId = int(doi)
        print 'searching doc of pubmed_id = %s' % pubmedId
        q = Docs.objects.filter(pubmed_id = pubmedId)

    except ValueError:
        print 'searching doc of doi = %s' % doi
        q = Docs.objects.filter(doi__exact = doi)

    if len(q):
        doc_id = q[0].pk
    else:
        print 'searchuin'
        q = Docs.objects.filter(pubmed_id = int(id))
        if len(q):
            doc_id = q[0].pk
        elif meta.get('doi'):
            q = Docs.objects.filter(doi__exact = meta['doi'])
            if len(q):
                doc_id = q[0].pk

    if doc_id:
        doc = q[0]
        journal = doc.journal
        arts = JournalArticles.objects.filter(pk=doc_id)
        art = None
        if len(arts):
            art = arts[0]
        if not meta['journal']['title']:
            meta['journal']['title'] = journal.title if journal else None
        if not meta['journal']['ISOAbbreviation']:
            meta['journal']['ISOAbbreviation'] = journal.iso_abbreviation if journal else None
        if not meta['journal']['issn']:
            meta['journal']['issn'] = journal.issn_print if journal else None
        if not meta['journal']['issn']:
            meta['journal']['issn'] = journal.issn_electronic if journal else None
        meta['journal']['volume'] = doc.volume
        meta['journal']['issue'] = doc.issue
        if not meta['journal']['pubDate']['year']:
            meta['journal']['pubDate']['year'] = art.year if art else None
        if not meta['journal']['pubDate']['month']:
            meta['journal']['pubDate']['month'] = art.month if art else None
        if not meta['journal']['pubDate']['day']:
            meta['journal']['pubDate']['day'] = art.day if art else None
        meta['journal']['pagination'] = art.pagination if art else None
        meta['first_page'] = doc.first_page
        meta['last_page'] = doc.last_page
        if not meta['title']:
            meta['title'] = doc.title
        if not meta['abstract']:
            meta['abstract'] = doc.abstract
        if not meta['authors']:
            meta['authors'] = doc.authors

    meta['doc_id'] = doc_id

    meta['chembl_like'] = "No"

    title =  urlquote(meta['title'])
    abstract  = urlquote(meta['abstract'])
    url = '%sCHEMBLLIKE/%s/%s' % (settings.PIPLINE_PILOT_ENDPOINT, title, abstract)
    try:
        result = requests.get(url, timeout=60)
        status = result.status_code

        if status != 200:
            pass
        else:
            if result.json()["Prediction"]:
                meta['chembl_like'] = "Yes"
    except:
        pass

    return meta

#-----------------------------------------------------------------------------------------------------------------------

def getStructure(mol):
    data = dict()
    if settings.OPEN_SOURCE:
        try:
            inchi = inchiFromPipe(mol, settings.INCHI_BINARIES_LOCATION['1.02'])
            data['InChI'] = inchi
            inchiKey = InchiToInchiKey(inchi)
            data['InChIKey'] = inchiKey
            smiles = smilesFromMol(mol)
            data['Canonical_Smiles'] = smiles
        except:
            pass
    else:
        url = '%scuration' % settings.PIPLINE_PILOT_ENDPOINT
        result = requests.post(url, data=mol, timeout=60)
        status = result.status_code

        if status != 200:
            raise Exception("URL %s has status %s for mol %s" % (url, status, mol))
        data = result.json()
    return data

#-----------------------------------------------------------------------------------------------------------------------

def InchiKeyFromMol(mol):
    return getStructure(mol).get('InChIKey', '')

#-----------------------------------------------------------------------------------------------------------------------

def getStatus():
    try:
        if settings.OPEN_SOURCE:
            return 6
        url = '%salive' % settings.PIPLINE_PILOT_ENDPOINT
        result = requests.get(url, timeout=60)
        status = result.status_code

        if status != 200:
            return 1
        res = result.json()
        if not res['up']:
            return 1   
        return int(res['resourceindex'])
    except:
        return 1    

#-----------------------------------------------------------------------------------------------------------------------

def ImageFromMolPP(mol, size):
    url = '%sctab2image/%s' % (settings.PIPLINE_PILOT_ENDPOINT, size)
    result = requests.post(url, data=mol, timeout=60)
    status = result.status_code

    if status != 200:
        raise Exception("URL %s has status %s for mol %s" % (url, status, mol))
    return result.json()['b64PNG']

#-----------------------------------------------------------------------------------------------------------------------

def ImageFromMol(mol, size=100):
    if settings.OPEN_SOURCE:
        molecule = Chem.MolFromMolBlock(str(mol))
        raw = Draw.MolToImage(molecule, size=(size, size))

        output = StringIO.StringIO()
        raw.save(output, 'PNG')
        ret = output.getvalue()
        output.close()
        return b64encode(ret)

    else:
        return ImageFromMolPP(mol, size)

#-----------------------------------------------------------------------------------------------------------------------

def cleanup(molfile, mode = 'cleanup'):
    #TODO: open source version of this function using RDKit should be implemented as well
    url = '%scleanup/%s' % (settings.PIPLINE_PILOT_ENDPOINT, mode)
    result = requests.post(url, data=molfile, timeout=60)
    status = result.status_code

    if status != 200:
        raise Exception("URL %s has status %s for mol %s" % (url, status, molfile))
    data = result.json()
    return data


#-----------------------------------------------------------------------------------------------------------------------

def compoundSearch(query, max_results = 5):
    import re
    from haystack.query import SearchQuerySet
    from chembl_business_model.models import MoleculeDictionary
    from chembl_business_model.models import ChemblIdLookup
    from chembl_business_model.models import CompoundStructures
    from chembl_business_model.models import Products
    from chembl_business_model.models import CompoundRecords
    from chembl_business_model.models import MoleculeSynonyms
    from chembl_business_model.models import Biotherapeutics

    if query.upper() == 'CHEMBL': # good luck with that
        return []

    if query.isdigit():
        try:
            q = MoleculeDictionary.objects.get(pk=int(query))
            return [{'label' : "%s (%s)" % (q.pref_name, q.chembl_id), 'value': int(query)}]
        except:
            pass

    if query.upper().startswith('CHEMBL'):
        try:
            q = ChemblIdLookup.objects.get(pk=query.upper())
            if q.entity_type == 'COMPOUND':
                mol = MoleculeDictionary.objects.get(chembl_id=query.upper())
                return [{'label' : "%s (%s)" % (mol.pref_name, query.upper()), 'value': mol.molregno}]
        except:
            q = SearchQuerySet().models(ChemblIdLookup).autocomplete(chembl_id=query)
            r = []
            if len(q):
                for chembl in q:
                    if chembl.object.entity_type == 'COMPOUND':
                        mol = MoleculeDictionary.objects.get(chembl_id=chembl.object.chembl_id)
                        r.append({'label' : "%s (%s)" % (mol.pref_name, query), 'value': mol.molregno})
                if len(r):
                    return r[:max_results]

    if len(query) == 27 and query[14] == '-' and query[25] == '-' and re.match('^([0-9A-Z\-]+)$',query):
        try:
            q = CompoundStructures.objects.get(standard_inchi_key=query).molecule
            return [{'label' : "%s (%s)" % (q.pref_name, q.chembl_id), 'value': q.molregno}]
        except:
            pass

    q = SearchQuerySet().models(MoleculeDictionary).autocomplete(pref_name=query)
    if len(q):
        r = []
        for mol in q:
            r.append({'label' : "%s (%s)" % (mol.object.pref_name, mol.object.chembl_id), 'value': mol.object.molregno})
        return r[:max_results]

    q = SearchQuerySet().models(Products).autocomplete(trade_name=query)
    if len(q):
        r = []
        for prod in q:
            mols = prod.object.moleculedictionary_set
            if mols.count():
                for mol in mols.all():
                    r.append({'label' : "%s (%s)" % (mol.pref_name, mol.chembl_id), 'value': mol.molregno})
        if len(r):
            return r[:max_results]

    q = SearchQuerySet().models(CompoundRecords).autocomplete(compound_name=query)
    if len(q):
        r = []
        for record in q:
            mol = record.object.molecule
            r.append({'label' : "%s (%s)" % (mol.pref_name, mol.chembl_id), 'value': mol.molregno})
        return r[:max_results]

    q = SearchQuerySet().models(CompoundRecords).autocomplete(compound_key=query)
    if len(q):
        r = []
        for record in q:
            mol = record.object.molecule
            r.append({'label' : "%s (%s)" % (mol.pref_name, mol.chembl_id), 'value': mol.molregno})
        return r[:max_results]

    q = SearchQuerySet().models(MoleculeSynonyms).autocomplete(synonyms=query)
    if len(q):
        r = []
        for synonym in q:
            mol = synonym.object.molecule
            r.append({'label' : "%s (%s)" % (mol.pref_name, mol.chembl_id), 'value': mol.molregno})
        return r[:max_results]

    q = SearchQuerySet().models(Biotherapeutics).autocomplete(description=query)
    if len(q):
        r = []
        for bio in q:
            mol = bio.object.molecule
            r.append({'label' : "%s (%s)" % (mol.pref_name, mol.chembl_id), 'value': mol.molregno})
        return r[:max_results]

#-----------------------------------------------------------------------------------------------------------------------

def md5Checksum(fh):
    fh.seek(0,0)
    m = hashlib.md5()
    while True:
        data = fh.read(8192)
        if not data:
            break
        m.update(data)
    return m.hexdigest()

#-----------------------------------------------------------------------------------------------------------------------
