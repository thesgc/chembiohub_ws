from django.conf import settings


from deployment.urls import urlpatterns as orignalpatterns

from cbh_chembl_ws_extension.urls import api_name
from tastypie.api import Api
from example_custom_deployment.resources import AlteredCompoundBatchResource

api = Api(api_name=api_name)

api.register(AlteredCompoundBatchResource())

print api.urls

urlpatterns = api.urls + orignalpatterns




