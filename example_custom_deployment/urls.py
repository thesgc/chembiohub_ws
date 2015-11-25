from django.conf.urls import patterns, include, url
from django.conf import settings

patt = ('',(r'^', include('cbh_chembl_ws_extension.urls')),
   (r'^', include('cbh_datastore_ws.urls')),
   
    (r'^django-rq/', include('django_rq.urls')),
)

if "tastypie_spore_docs" in settings.INSTALLED_APPS:
    patt += ((r'^', include('tastypie_spore_docs.urls'))),


#Here we generate some override URLs

DEFAULT_API_NAME = 'chemblws'

try:
    api_name = settings.WEBSERVICES_NAME
except AttributeError:
    api_name = DEFAULT_API_NAME


from cbh_chembl_ws_extension.urls import api


urlpatterns = patterns(

    *patt
    )



