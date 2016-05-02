"""
Background task functions which are called from elsewhere using function signatures and async_iter
Tasks are run using django_q
"""

from cbh_chem_api.resources import index_batches_in_new_index
from cbh_chembl_model_extension.models import generate_structure_and_dictionary, CBHCompoundBatch
from chembl_business_model.models import CompoundMols
from django.contrib.postgres.aggregates import ArrayAgg
from cbh_core_model.models import Project




def process_batch_list(batch_list):
    """Generate the structure and dictionary for a list of compound batches"""
    batches = [generate_structure_and_dictionary(batch) for batch in batch_list]
    index_batches_in_new_index(batches)
    return batches




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



