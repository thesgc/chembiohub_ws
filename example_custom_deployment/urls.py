from django.conf import settings


from deployment.urls_v2 import urlpatterns as orignalpatterns

from cbh_chem_api.urls import api_name
from tastypie.api import Api
from example_custom_deployment.resources import AlteredCompoundBatchResource

api = Api(api_name=api_name)

api.register(AlteredCompoundBatchResource())

print api.urls

urlpatterns = api.urls + orignalpatterns




