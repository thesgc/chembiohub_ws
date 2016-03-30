from cbh_chem_api.resources import index_batches_in_new_index
from cbh_chembl_model_extension.models import generate_structure_and_dictionary, CBHCompoundBatch
from chembl_business_model.models import CompoundMols

def process_batch_list(batch_list):
    batches = [generate_structure_and_dictionary(batch) for batch in batch_list]
    index_batches_in_new_index(batches)
    return batches

def get_structure_search_for_project(project_id, search_type, smiles):
    queryset = CompoundMols.objects.filter(molecule__project_id=project_id)
    search_function = getattr(queryset, search_type)
    batch_ids = search_function(smiles).values_list("molecule_id", flat=True)
    batch_ids = CBHCompoundBatch.objects.filter(related_molregno_id__in=batch_ids).values_list("id", flat=True)
    return {"project_id" : project_id,"batch_ids" : batch_ids}