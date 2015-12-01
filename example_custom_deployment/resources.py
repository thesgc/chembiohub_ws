from cbh_chembl_ws_extension.compounds import CBHCompoundBatchResource
from django.http import HttpRequest

class AlteredCompoundBatchResource(CBHCompoundBatchResource):
    def validate_multi_batch(self, multiple_batch, bundle, request, batches):
        response = super(AlteredCompoundBatchResource, self).validate_multi_batch(multiple_batch, bundle, request, batches)
        
        request.GET = request.GET.copy()
        request.GET["format"] = "sdf"
        request.GET["current_batch"] = multiple_batch.id
        sdf_file_resp = super(AlteredCompoundBatchResource, self).get_part_processed_multiple_batch(request)
        print sdf_file_resp.__dict__
        return response
