from django.conf.urls import patterns, url
from views import UploadView


# JSON REQUESTS
urlpatterns = patterns('',
    url(r'^upload/(?P<project_id>[0-9]+)/$', UploadView.as_view(), name="flowv2_upload"),
)
