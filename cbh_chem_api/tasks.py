"""
Background task functions which are called from elsewhere using function signatures and async_iter
Tasks are run using django_q
"""

import re
from cbh_chembl_model_extension.models import generate_structure_and_dictionary, CBHCompoundBatch
from chembl_business_model.models import CompoundMols
from django.contrib.postgres.aggregates import ArrayAgg
from cbh_core_model.models import Project
from rdkit import Chem
from rdkit.Chem.AllChem import Compute2DCoords
from cbh_chem_api import chemdraw_reaction
from cbh_utils.parser import parse_pandas_record, parse_sdf_record, apply_json_patch, get_uncurated_fields_from_file, get_all_sdf_headers
import itertools
import math
import pandas as pd
import copy
from difflib import SequenceMatcher as fuzzymatch
import shortuuid
from django.conf import settings
from cbh_chembl_model_extension.models import MoleculeDictionary
from cbh_utils import elasticsearch_client
from pybel import readfile, readstring
from django_q.tasks import async_iter, result, async
from cbh_chem_api.resources import index_batches_in_new_index


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

def process_file_request(
                         multiple_batch,
                         bundledata, 
                         creator_user,
                         schemaform,
                         correct_file,
                         session_key):
    from cbh_chem_api.compounds  import CBHCompoundUploadResource
    cbr_instance = CBHCompoundUploadResource()
    batches = []
    headers = []
    errors = []
    fielderrors = {}
    fielderrors["stringdate"] = set([])
    fielderrors["number"] = set([])
    fielderrors["integer"] = set([])
    structure_col = bundledata.get("struccol", "")
    
    if (".cdx" in correct_file.extension):
        mols = [mol for mol in readfile(
            str(correct_file.extension[1:]), str(correct_file.file.name), )]
        rxn = None
        index = 0
        if correct_file.extension == '.cdxml':
            # Look for a stoichiometry table in the reaction file
            rxn = chemdraw_reaction.parse(str(correct_file.file.name))
            headers = ["%Completion",
                       "%Yield",
                       "Expected Moles",
                       "Product Moles",
                       "Expected Mass",
                       "Product Mass",
                       "MW",
                       "role",
                       "Purity",
                       "Limit Moles",
                       # "Formula",
                       "Equivalents",
                       "Measured Mass"]
        multiple_batch.batch_count = len(mols)
        for pybelmol in mols:
            molfile = pybelmol.write("mdl")
            if molfile.strip() and molfile != "*":
                rd_mol = Chem.MolFromMolBlock(molfile, sanitize=False)
                '''
                Pybel can read some hypervalent molecules that RDKit cannot read
                Therefore currently these molecules are outputted as images and sent back to the front end
                https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg04466.html
                '''
                if rd_mol:
                    smiles = Chem.MolToSmiles(rd_mol)
                    if smiles.strip():
                        try:
                            b = CBHCompoundBatch.objects.from_rd_mol(
                                rd_mol, orig_ctab=molfile, smiles=smiles, project=bundledata["project"])
                        except Exception, e:
                            b = None
                            index = index - 1
                        if b:
                            if rxn:
                                # Here we set the uncurated fields equal to
                                # the reaction data extracted from Chemdraw
                                b.uncurated_fields = rxn.get(
                                    pybelmol.title, {})
                            batches.append(b)
                        else:
                            errors.append({"index": index+1, "image": pybelmol.write(
                                "svg"), "message": "Unable to produce inchi from this molecule"})
                else:
                    errors.append({"index": index+1, "image": pybelmol.write(
                        "svg"), "message": "Invalid valency or other error parsing this molecule"})
                index += 1

    else:
        if (correct_file.extension == ".sdf"):
            headers = get_all_sdf_headers(correct_file.path)
            uncurated, ctabs, ctab_parts = get_uncurated_fields_from_file(correct_file, fielderrors)
            multiple_batch.batch_count = len(ctabs)
            multiple_batch.save()
            args = [(index, arguments[0], arguments[1], arguments[2], bundledata["project"]) for index, arguments in enumerate(itertools.izip( ctabs, ctab_parts,  uncurated))]
            #split data into 3 parts
            list_size = math.ceil(float(len(args))/2.0)

            arg_chunks = [c for c in chunks(args, list_size)]
            id = async_iter("cbh_chem_api.tasks.get_batch_from_sdf_chunks", arg_chunks)
            lists_of_batches =  result(id, wait=1000000)
            for batchlist in lists_of_batches:
                if isinstance(batchlist, basestring):
                    
                    raise Exception(batchlist)
            batches = [inner for outer in lists_of_batches for inner in outer]


            # batches = get_batch_from_sdf_chunks(args) 
                
        elif(correct_file.extension in (".xls", ".xlsx")):
            # we need to know which column contains structural info - this needs to be defined on the mapping page and passed here
            # read in the specified structural column

            headerswithdata = set([])

            df = None
            df = pd.read_excel(correct_file.file)
            

            # read the smiles string value out of here, when we know which
            # column it is.
            

            row_iterator = df.iterrows()
            headers = list(df)
            headers = [h.replace(".", "__") for h in headers]
            headers = [h.replace("/", "__") for h in headers]
            df.columns = headers
            # Only automap on the first attempt at mapping the smiles
                # column
            if not structure_col and not bundledata.get("headers", None):
                max_score = 0
                for header in headers:
                    # fuzzy matching for smiles - this should also
                    # match things like "canonical_smiles"
                    hdr = re.sub('[^0-9a-zA-Z]+', ' ', header)
                    for h in hdr.split(" "):
                        h = h.strip()
                        if h:
                            score = fuzzymatch(
                                a="smiles", b=h.lower()).ratio()
                            if score > max_score and score > 0.9:
                                structure_col = header
                                max_score = score
                                automapped_structure = True
            
            args = [(index, row, structure_col, bundledata["project"]) for index, row in row_iterator]

            list_size = math.ceil(float(len(args))/2.0)

            arg_chunks = [c for c in chunks(args, list_size)]
            id = async_iter("cbh_chem_api.tasks.get_batch_from_xls_chunks", arg_chunks)
            lists_of_batches =  result(id, wait=1000000)
            for batchlist in lists_of_batches:
                if isinstance(batchlist, basestring):
                    raise Exception(batchlist)
            batches = [inner for outer in lists_of_batches for inner in outer]
          

            for b, args in itertools.izip(batches, args):
                if dict(b.uncurated_fields) == {}:
                # Only rebuild the uncurated fields if this has not
                # been done before
                    parse_pandas_record(
                        headers, b, "uncurated_fields", args[1], fielderrors, headerswithdata)
               
            headers = [hdr for hdr in headers if hdr in headerswithdata]
            
    
    for b in batches:
        if b:
            b.multiple_batch_id = multiple_batch.pk
            b.created_by = creator_user.username
            b.created_by_id = creator_user.id

    bundledata["fileerrors"] = errors
    bundledata["automapped"] = 0
    
    if not bundledata.get("headers", None):
        bundledata["headers"] = []
        for header in headers:
            copyto = ""
            automapped = False
            operations = []
            if header == structure_col:
                copyto = "SMILES for chemical structures"
                if automapped_structure:
                    automapped = True
            else:
                form = copy.deepcopy(schemaform)
                copyto = ""
                max_score = 0
                for form_item in form:
                    score = fuzzymatch(
                        a=form_item["key"].lower(), b=header.lower()).ratio()
                    if score > max_score and score > 0.9:
                        matched_item = form_item
                        copyto = matched_item["key"]
                        automapped = True
                        if(matched_item.get("field_type", "") == "uiselecttags"):
                            operations.append(
                                {"op": "split", "path": "/uncurated_fields/" + header})
                            operations.append(
                                {"op": "move", "path": "/custom_fields/" + matched_item["key"], "from": "/uncurated_fields/" + header})
                        else:
                            operation = {
                                "op": "move", "path": "/custom_fields/" + matched_item["key"], "from": "/uncurated_fields/" + header}
                            operations.append(operation)
                            if(matched_item.get("format", "") == "date"):
                                operations.append(
                                    {"op": "convertdate", "path": "/custom_fields/" + matched_item["key"]})
                        #set the max score so less well matched content than this is ignored
                        max_score = score

            bundledata["headers"].append({
                "name": header,
                "automapped": automapped,
                "copyto": copyto,
                "operations": operations,
                "fieldErrors": {
                    "stringdate": header in fielderrors["stringdate"],
                    "integer": header in fielderrors["integer"],
                    "number": header in fielderrors["number"]
                }
            })

    validate_multi_batch(cbr_instance, multiple_batch, bundledata, session_key, batches)
    return True

    
def validate_multi_batch(cbr_instance, multiple_batch, bundledata, session_key, batches):
    """Generate a set of staticstics about a set of data that has been uploaded"""
    batches_not_errors = [batch for batch in batches if batch and not batch.warnings.get(
        "parseerror", None) and not batch.warnings.get("smilesParseError", None)]


    for b in batches_not_errors:
        b.properties["action"] = "New Batch"
    batches_with_structures = [
        batch for batch in batches_not_errors if batch.ctab]
    blinded_data = [
        batch for batch in batches_not_errors if not batch.ctab]
    sdfstrings = [batch.ctab for batch in batches_with_structures]
    sdf = "\n".join(sdfstrings)

    filename = "/tmp/" + shortuuid.ShortUUID().random()
    text_file = open(filename, "w")
    text_file.write(sdf)
    text_file.close()
    from subprocess import PIPE, Popen
    p = Popen([settings.INCHI_BINARIES_LOCATION['1.02'],
               "-STDIO",  filename], stdout=PIPE, stderr=PIPE)

    a = p.communicate()
    inchis = {}

    #PB - there is an assumption here that everything that has a structure will generate an inChi without issue. This is not the case.
    #Where a molecule does not generate an inchi, there will be a key error looking up the inchi in inchiparts, as anything that cannot 
    #generate an inchi will be missing from inchiparts, i.e. 50 structures with 1 error will have 49 entries in inchiparts, and this 
    #will in turn bin the whole file - not great when we can handle erroring structures elsewhere

    error_locs = []

    #a[0] holds the generated inchis. a[1] holds all of the error and warning information (if any)
    errorparts = a[1].split("\nError")
    if(len(errorparts) > 1):
        for i, errorp in enumerate(errorparts):
            #split on 'structure #', then get the number given
            if(i > 0):
                splits = errorp.split('structure #')

                error_loc = splits[1].split('.')[0]
                #convert to number, put this number in an errors list
                error_locs.append(error_loc)

    err_batches = []
    #for the errors found, remove from non-error lists and flag as erroring
    for error_no in error_locs:
        error_no_int = int(float(error_no)) - 1

        #find structures at the position indicated - 1 (for 0-indexed list)
        err_batch = batches_with_structures[error_no_int]
        err_batches.append(err_batch)

    #we can't remove these while looping through err_locs as it messes up the list order and gives arrayindex exceptions
    for err_batch in err_batches:

        #remove from batches_with_structures and batches_not_errors
        batches_with_structures.remove(err_batch)
        batches_not_errors.remove(err_batch)

        #flag this batch as erroring due to inability to generate anything for the standard_inchi_key field
        batches_index = batches.index(err_batch)
        batches[batches_index].warnings["inchicreationerror"] = "true"
        batches[batches_index].properties["action"] = "Ignore"


    inchiparts = a[0].split("\nStructure:")

    for i, inch in enumerate(inchiparts):
        parts = inch.split("\n")
        if len(parts) == 1:
            continue
        ints = [s for s in parts[0].split() if s.isdigit()]
        part = "".join(ints)
        inchis[part] = parts[1]
    if not bundledata.get("fileerrors"):
        bundledata["fileerrors"] = []
    new_uploaded_data = []
    already_found = set([])
    duplicates = set([])
    for i, batch in enumerate(batches_with_structures):
        if (str(i+1) in error_locs):
            batch.standard_inchi = None
        else: 
            batch.standard_inchi = inchis[str(i+1)]
        batch.validate(temp_props=False)
        if batch.standard_inchi_key in already_found:
            # setting this in case we change it later
            duplicates.add(batch.standard_inchi_key)
        else:
            already_found.add(batch.standard_inchi_key)

        new_uploaded_data.append(batch)
    already_in_db = MoleculeDictionary.objects.filter(project=bundledata[
                                                      "project"], structure_type="MOL", structure_key__in=already_found).values_list("structure_key", flat=True)
    already_in_db = set(already_in_db)

    bundledata["new"] = 0
    new_data = set([])
    duplicate_overlaps = set([])
    duplicate_new = set([])
    for batch in batches_with_structures:
        if batch.standard_inchi_key in duplicates:
            batch.warnings["duplicate"] = True
        if batch.standard_inchi_key in already_in_db:
            batch.warnings["overlap"] = True
            if batch.standard_inchi_key in duplicates:
                batch.warnings["duplicate"] = True
                duplicate_overlaps.add(batch.standard_inchi_key)
        else:
            batch.warnings["new"] = True

            new_data.add(batch.standard_inchi_key)
            if batch.standard_inchi_key in duplicates:
                batch.warnings["duplicate"] = True
                duplicate_new.add(batch.standard_inchi_key)

    for batch in batches_with_structures:
        if batch.warnings.get("withoutstructure") == True:
            del batch.warnings["withoutstructure"]
    for batch in blinded_data:
        batch.warnings["withoutstructure"] = True
    bundledata["batchstats"] = {}
    bundledata["batchstats"]["withstructure"] = len(
        batches_with_structures)
    bundledata["batchstats"]["parseerrors"] = len(batches) - len(batches_not_errors) + len(
        [b for b in batches_not_errors if b.warnings.get("parseerror", False) == "true"])
    bundledata["batchstats"]["withoutstructure"] = len(blinded_data)
    bundledata["batchstats"]["total"] = len(batches)
    bundledata["compoundstats"] = {}
    bundledata["compoundstats"]["total"] = len(
        already_in_db) + len(new_data)
    bundledata["compoundstats"]["overlaps"] = len(already_in_db)
    bundledata["compoundstats"]["new"] = len(new_data)
    bundledata["compoundstats"][
        "duplicateoverlaps"] = len(duplicate_overlaps)
    bundledata["compoundstats"]["duplicatenew"] = len(duplicate_new)
    bundledata["multiplebatch"] = multiple_batch.pk


    cbr_instance.set_cached_temporary_batches(
        batches, multiple_batch.id, session_key)
    
    #bundledata["objects"] = fifty_batches_for_first_page
    index_name = elasticsearch_client.get_temp_index_name(
        session_key, multiple_batch.id)
    elasticsearch_client.get_action_totals(index_name, bundledata)
    multiple_batch.uploaded_data = bundledata
    multiple_batch.save()





def get_batch_from_sdf_chunks(lists):
    batches = []
    for args in lists:
        batches.append(get_batch_from_sdf(*args))
    return batches

def get_batch_from_sdf(index, sdf_data, ctab_part, uncur, project):
    mol = Chem.MolFromMolBlock(ctab_part)
    if mol is None:
        b = CBHCompoundBatch.objects.blinded(
            project=project)
        b.warnings["parseerror"] = "true"
        b.properties["action"] = "Ignore"

        b.uncurated_fields = uncur

    else:
        try:
            b = CBHCompoundBatch.objects.from_rd_mol(
                mol, orig_ctab=sdf_data, project=project)
        except Exception, e:
            b = CBHCompoundBatch.objects.blinded(
            project=project)
            b.warnings["parseerror"] = "true"
            b.properties["action"] = "Ignore"

        b.uncurated_fields = uncur
    return b



def get_batch_from_xls_chunks(args_list):
    batches = []
    for args in args_list:
        index, row, structure_col, project = args
        
        batch = get_batch_from_xls_row(index, row, structure_col, project)
        batches.append(batch)
    return batches

def get_batch_from_xls_row(index, row, structure_col, project):
    if structure_col:
        smiles_str = row[structure_col]
        try:
            struc = Chem.MolFromSmiles(smiles_str)
            if struc:
                Compute2DCoords(struc)
                b = CBHCompoundBatch.objects.from_rd_mol(
                    struc, smiles=smiles_str, project=project, reDraw=True)
                b.blinded_batch_id = None
            else:
                raise Exception("Smiles not processed")

        except Exception, e:
            b = CBHCompoundBatch.objects.blinded(
                project=project)
            b.original_smiles = smiles_str
            b.warnings["parseerror"] = "true"
            b.properties["action"] = "Ignore"
    else:
        b = CBHCompoundBatch.objects.blinded(
            project=project)
    return b





def get_structure_search_for_projects(project_ids, search_type, smiles):
    """Run a substructure or exact match search across the requested projects via the manager
    methods provided on the compoundstructures method in chembl_core_db"""
    queryset = CompoundMols.objects.filter(molecule__project_id__in=project_ids)
    search_function = getattr(queryset, search_type)
    mol_ids = search_function(smiles).order_by("-molecule_id").values_list("molecule_id", flat=True)

    Project.objects.filter(cbhcompoundbatch__related_molregno_id__in=mol_ids)

    batch_ids = CBHCompoundBatch.objects.filter(project_id__in=project_ids)\
                                        .filter(related_molregno_id__in=mol_ids)\
                                        .order_by("project_id")\
                                        .values("project_id")\
                                        .annotate(batch_ids=ArrayAgg("id"))


    return batch_ids



def save_multiple_batch( 
                        multiple_batch, 
                        creator_user, 
                        session_key):

    from cbh_chem_api.compounds  import CBHCompoundUploadResource
    cbr_instance = CBHCompoundUploadResource()

    limit = math.ceil(float(multiple_batch.batch_count)/2.0)
    offset = 0
    batches = []
    hasMoreData = True
    

    datasets = []
    for run in range(0,3):
        datasets.append(( 
                        multiple_batch, 
                        creator_user, 
                        session_key,
                        limit, 
                        offset))
        offset += limit
    id = async_iter("cbh_chem_api.tasks.process_batch_list", datasets)
    lists_of_batches = result(id, wait=1000000)
    batches = [inner for outer in lists_of_batches for inner in outer]
    if multiple_batch.uploaded_file:
        cbr_instance.alter_batch_data_after_save( batches , multiple_batch.uploaded_file.file , multiple_batch)
    index_batches_in_new_index(batches)
    elasticsearch_client.delete_index(
        elasticsearch_client.get_temp_index_name(session_key, multiple_batch.id))
    print multiple_batch.project_id
    #cbr_instance.after_save_and_index_hook(request, id, multiple_batch.project_id)
    return True





def process_batch_list(multiple_batch, 
                        creator_user, 
                        session_key,
                        limit, 
                        offset):
    """Generate the structure and dictionary for a list of compound batches"""
    from cbh_chem_api.compounds  import CBHCompoundUploadResource
    cbr_instance = CBHCompoundUploadResource() 


    bundles = cbr_instance.get_cached_temporary_batch_data(
        multiple_batch.id,  {"query": '{"term" : {"properties.action.raw" : "new batch"}}',"limit": limit, "offset": offset, "sorts": '[{"id": {"order": "asc"}}]'}, session_key)
   # allowing setting of headers to be fale during saving of drawn molecule
    if multiple_batch.uploaded_data["headers"]:
        for d in bundles["objects"]:
            cbr_instance.patch_dict(d, copy.deepcopy(multiple_batch.uploaded_data["headers"]))
    set_of_batches = cbr_instance.get_cached_temporary_batches(
        bundles, None,  bundledata=multiple_batch.uploaded_data)
    for batch in set_of_batches["objects"]:
        if(batch.obj.properties.get("action", "") == "New Batch"):
            batch.obj.created_by = creator_user.username
            batch.obj.created_by_id = creator_user.id
            batch.obj.id = None
            batch.obj.multiple_batch_id = multiple_batch.id

            # batch.multi_batch_id = id
    batches = [b.obj for b in set_of_batches["objects"]]

    

    batches = [generate_structure_and_dictionary(batch) for batch in batches]
    
    return batches


