from cbh_chem_api.resources import index_batches_in_new_index
from cbh_chembl_model_extension.models import generate_structure_and_dictionary


def process_batch_list(batch_list):
    batches = [generate_structure_and_dictionary(batch) for batch in batch_list]
    index_batches_in_new_index(batches)
    return batches

