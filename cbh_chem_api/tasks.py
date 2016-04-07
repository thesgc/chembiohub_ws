from cbh_chem_api.resources import index_batches_in_new_index
from cbh_chembl_model_extension.models import generate_structure_and_dictionary, CBHCompoundBatch
from chembl_business_model.models import CompoundMols




def process_batch_list(batch_list):
    batches = [generate_structure_and_dictionary(batch) for batch in batch_list]
    index_batches_in_new_index(batches)
    return batches

def is_simple(smiles):
    if len(smiles) < 5:
        return True
    return False


def get_structure_search_for_projects(project_ids, search_type, smiles):
    queryset = CompoundMols.objects.filter(molecule__project_id__in=project_ids)
    search_function = getattr(queryset, search_type)
    mol_ids = search_function(smiles).order_by("-molecule_id").values_list("molecule_id", flat=True)
    restricted = False
    if search_type == "with_substructure" and is_simple(smiles):

        mol_ids = mol_ids[:1000]
        mol_ids = list(mol_ids)
        if len(mol_ids) == 1000:
            restricted = True

    batch_ids = CBHCompoundBatch.objects.filter(related_molregno_id__in=mol_ids).values_list("id", flat=True)

    return {"project_ids" : project_ids,"batch_ids" : batch_ids, "restricted" : restricted}



