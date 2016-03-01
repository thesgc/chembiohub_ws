from cbh_core_api.authorization import ProjectAuthorization
from cbh_core_api.elasticsearch_client import get_project_compound_index_name

class CompoundBatchAuthorization(ProjectAuthorization):
    def get_list_of_authorized_indices(self, request):
        restricted_and_unrestricted_projects = get_projects_where_fields_restricted(request.user)
        return [elasticsearch_client.get_project_compound_index_name(project_id, restriction) 
                    for project_id, restriction 
                    in restricted_and_unrestricted_projects.items()]


