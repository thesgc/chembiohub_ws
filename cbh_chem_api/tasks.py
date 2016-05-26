"""
Background task functions which are called from elsewhere using function signatures and async_iter
Tasks are run using django_q
"""


from cbh_chembl_model_extension.models import generate_structure_and_dictionary, CBHCompoundBatch
from chembl_business_model.models import CompoundMols
from django.contrib.postgres.aggregates import ArrayAgg
from cbh_core_model.models import Project
from rdkit import Chem
from rdkit.Chem.AllChem import Compute2DCoords



def process_batch_list(batch_list):
    """Generate the structure and dictionary for a list of compound batches"""
    batches = [generate_structure_and_dictionary(batch) for batch in batch_list]
    
    return batches



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



