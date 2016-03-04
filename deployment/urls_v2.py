from django.conf.urls import patterns, include, url
from django.conf import settings

patt = ('',(r'^', include('cbh_chem_api.urls')),
   
    (r'^django-rq/', include('django_rq.urls')),
    
)

if "tastypie_spore_docs" in settings.INSTALLED_APPS:
    patt += ((r'^', include('tastypie_spore_docs.urls'))),

urlpatterns = patterns(

    *patt
    )
