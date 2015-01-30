from cbh_chembl_ws_extension import __version__ as ws_version
from django.conf.urls import patterns, include, url
from django.conf import settings

urlpatterns = patterns('',

   (r'^', include('cbh_chembl_ws_extension.urls')),

    )

