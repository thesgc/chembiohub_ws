from django.conf.urls import patterns, include, url
from django.conf import settings

urlpatterns = patterns('',

   (r'^', include('cbh_chembl_ws_extension.urls')),
   (r'^', include('cbh_datastore_ws.urls')),
   (r'^', include('tastypie_spore_docs.urls')),

    )

