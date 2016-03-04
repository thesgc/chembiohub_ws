from django import http
from django import forms
from django.conf import settings
from django.views.generic.base import View
from django.shortcuts import get_object_or_404
from cbh_core_model.models import CBHFlowFile as FlowFile, CBHFlowFileChunk as FlowFileChunk
from tastypie.resources import ModelResource, Resource, ALL, ALL_WITH_RELATIONS
#from tastypie.authorization import Authorization
from tastypie.authentication import SessionAuthentication
from tastypie.serializers import Serializer
import mimetypes


class FlowFileForm(forms.Form):
    file = forms.FileField()


class UploadView(View):
    def dispatch(self, request, *args, **kwargs):
        # get flow variables
        self.flowChunkNumber = int(request.REQUEST.get('flowChunkNumber'))
        self.flowChunckSize = int(request.REQUEST.get('flowChunkSize'))
        self.flowCurrentChunkSize = int(request.REQUEST.get('flowCurrentChunkSize'))
        self.flowTotalSize = int(request.REQUEST.get('flowTotalSize'))
        self.flowIdentifier = request.REQUEST.get('flowIdentifier')
        self.flowFilename = request.REQUEST.get('flowFilename')
        self.flowRelativePath = request.REQUEST.get('flowRelativePath')
        self.flowTotalChunks = int(request.REQUEST.get('flowTotalChunks'))


        # identifier is a combination of session key and flow identifier
        # create a new identifier which uses the project ID also
        print(request)
        self.identifier = ('%s-%s' % (request.session.session_key, self.flowIdentifier))[:255]
        return super(UploadView, self).dispatch(request, *args, **kwargs)

    def get(self, *args, **kwargs):
        """
        Flow.js test if chunk exist before upload it again.
        Return 200 if exist.
        """
        get_object_or_404(FlowFileChunk, number=self.flowChunkNumber, parent__identifier=self.identifier)
        return http.HttpResponse(self.identifier)

    def post(self, request, *args, **kwargs):
        """
        Upload the file by chunks
        """

        # get file or create if doesn't exist the identifier
        flow_file, created = FlowFile.objects.get_or_create(identifier=self.identifier, defaults={
            'original_filename': self.flowFilename,
            'total_size': self.flowTotalSize,
            'total_chunks': self.flowTotalChunks,
            'project_id': kwargs.get('project_id'),
        })

        # validate the file form
        form = FlowFileForm(request.POST, request.FILES)
        if not form.is_valid():
            return http.HttpResponseBadRequest(form.errors)

        # avoiding duplicated chucks
        chunk, created = flow_file.chunks.get_or_create(number=self.flowChunkNumber, defaults={
            'file': form.cleaned_data['file'],
        })
        if not created:
            chunk.file = form.file
            chunk.size = form.size
            chunk.save()

        return http.HttpResponse(flow_file.identifier)

#this is for requesting identifiers
class CBHFlowFileResource(ModelResource):
    #sheet_names = fields.ListField()

    class Meta:
        detail_uri_name = 'identifier'
        # Must be false to not give secret key away
        include_resource_uri = True
        allowed_methods = ['get', ]
        resource_name = 'cbh_flowfiles'
        queryset = FlowFile.objects.all()
        #filtering = {"identifier": ALL_WITH_RELATIONS}
        exclude = ['identifier']

    # def dehydrate_sheet_names(self, bundle):
    #     try:
    #         if bundle.obj.extension == "xlsx":
    #             sheetnames = get_sheetnames(bundle.obj.path)
    #             return sheetnames
    #         return []
    #     except XLRDError:
    #         raise BadRequest("Incorrect format or corrupt file")
    def alter_detail_data_to_serialize(self, request, data):
        sessionid = data.request.COOKIES.get(settings.SESSION_COOKIE_NAME, "None")
        print data.data
        bits = data.data["resource_uri"].split("/")
        bits[3] = bits[3][len(sessionid) + 1:]
        data.data["resource_uri"] = "/".join(bits)
        del data.data["identifier"]
        downloadres = CBHFlowFileDownloadResource()
        data.data["download_uri"] = downloadres._build_reverse_url('api_dispatch_detail', kwargs={"api_name": downloadres._meta.api_name, "resource_name": downloadres._meta.resource_name,  "pk": data.data["id"]})
        return data

    def get_list(self, *args, **kwargs):
        raise NotImplemented

    def obj_get(self, bundle, **applicable_filters):
        """
        An ORM-specific implementation of ``apply_filters``.
        The default simply applies the ``applicable_filters`` as ``**kwargs``,
        but should make it possible to do more advanced things.
        """
        if applicable_filters.get("identifier", None):
            print("applicable filters not none")
            applicable_filters["identifier"] = "%s-%s" % (
                bundle.request.COOKIES.get(settings.SESSION_COOKIE_NAME, "None"), applicable_filters["identifier"])
        return super(CBHFlowFileResource, self).obj_get(bundle, **applicable_filters)

#this is used to make the download request
class CBHFlowFileDownloadResource(ModelResource):

    class Meta:
        always_return_data = True  # required to add the elasticsearch data
        resource_name = 'cbh_downloads'
        default_format = 'application/json'
        include_resource_uri = True
        allowed_methods = ['post', 'get']
        default_format = 'application/json'
        serializer = Serializer()
        authentication = SessionAuthentication()
        queryset = FlowFile.objects.all()
        #need to replace this authorization with a custom one to check user can access the DPC
        #authorization = AttachmentAuthorization()
    #this is to render
    def get_detail(self, request, **kwargs):
        """
        Returns a single serialized resource.
        Calls ``cached_obj_get/obj_get`` to provide the data, then handles that result
        set and serializes it.
        Should return a HttpResponse (200 OK).
        """
        basic_bundle = self.build_bundle(request=request)

        obj = self.obj_get(basic_bundle, **kwargs)

        bundle = self.build_bundle(obj=obj, request=request)
        bundle = self.full_dehydrate(bundle)
        bundle = self.alter_detail_data_to_serialize(request, bundle)

        #return our response here
        #get extension from the FlowFile object
        #match this to a dictionary of mimetypes with extensions
        fb = obj.file.read()
        mimetype = mimetypes.guess_type(obj.full_path)[0]
        print(obj.original_filename)
        #if mimetype.index('spreadsheetml') > 0:
        
        response = http.HttpResponse(fb, content_type=mimetypes.guess_type(obj.full_path)[0])
        #if it's not an image, it is a download link - add the necessary content disposition info
        if(mimetype.count('image') == 0):
            response['Content-Disposition'] = 'attachment; filename=%s' % obj.original_filename
        return response
