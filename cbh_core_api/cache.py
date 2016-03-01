from django.core.cache import cache
from django.http import HttpResponse
from django.conf import settings



class CachedResource(object):

    """
    CachedResource Mixin

    This mixin allows TastyPie to cache the response at the
    begining of the request cycle.  It attempts to find a cache key
    that matches the API call being made.  If found it will return
    that cached response instead of making the full API rounds.

    This uses the standard cache defined in settings.py

    Most methods need overriding, but we'll try to keep things
    simple where we can with defaults.

    ---- General Flow Explained ----

    * Override the dispatch method of TastyPie. This is where most
    of the magic happens.
    * After making some sanity checks we'll look for a cache
    * To look for the cache we must get a key from the request params
    and compare it to the existing cache to see if one exists.
    * We must use the same format to return the cache because we are
    cacheing the entire serialized response.
    * If a cache is found, return it,  if not go through the request
    cycle and eventually store a cache using the CachedSerializer() class.
    * The Cached serializer class just makes caches and does not other
    altering of the request/response cycle.

    ---- End ----


    """

    valid_cache_formats = ['json', 'jsonp', 'xml', ]
    default_params = {'limit': 20}
    content_type_dict = {'json': 'application/json', 'jsonp':
                         'application/javascript', 'xml': 'application/xml'}
    valid_cache_get_keys = ['format', 'limit']

    def _prepare_params(self, params):
        """
        The serialzer expects to see some defaults for
         things like Limit.
        """
        for k, v in self.default_params.iteritems():
            if k not in params:
                params[k] = ['%s' % self.default_params[k]]
        return params

    def is_jsonp(self, request):
        """
        Determines if the request is of the JSONP format.  We'll
        need this later in order to determine if we need a callback
        in the response.
        """
        if request.GET.get('format') == 'jsonp':
            return True
        return False

    def make_cache_key(self, params):
        """
        Takes a Dictionary of query parameters and returns
        a String key.

        TODO: Is the resource specific or can this be global for all pulse apps?
        """

        key = ''

        for param_key in self.valid_cache_get_keys:
            # We need to check to make sure the order is the same
            # for both set and get keys.
            value = params.get(param_key, None)
            if value:
                key += '%s-%s-' % (param_key, value[0])

        return key

    def get_cache_key(self, params):
        """
        Make a cache key given a dict.  Can either call one
        from a util module or define it here.
        """
        return self.make_cache_key(params)

    def _set_cache(self, key, data):
        """ Set the Pulse API Front end Cache. """
        cache.set(key, data, getattr(settings, 'API_CACHE_LENGTH', 900))

    def get_cache(self, request, force_jsonp=False):
        """
        Retrieve a cache of the response for this query
        or none if no cache is found.
        """
        if self.valid_cache_format(request):
            params = dict(request.GET.iterlists())
            if force_jsonp:
                params['format'] = ['json']
            key = self.get_cache_key(self._prepare_params(params))
            data = cache.get(key)
            return data
        return False

    def valid_cache_format(self, request):
        """ Checks to see if the request
        is a format that we have cache support
        for.
        """
        return request.GET.get('format', None) not in ["xlsx", "sdf", "csv"]

    def _get_valid_content_type(self, format):
        """ Return the proper content type for the given
        format.  Uses the content_type_dict.
        """
        return self.content_type_dict[format]

    def create_response(self, request, data, response_class=HttpResponse, **response_kwargs):
        """
        Override of the default create_response method in resources.py.  This grabs the serialized response
        and caches it before passing back to the default method.
        """
        params = dict(request.GET.iterlists())
        if self.is_jsonp(request):
            params['format'] = ['json']
        key = self.get_cache_key(self._prepare_params(params))
        desired_format = self.determine_format(request)
        if self.is_jsonp(request):
            # if this is JSONP follow through if it was JSON
            # then we will wrap it with a call back later.
            desired_format = self.content_type_dict['json']
        self._set_cache(key, self.serialize(request, data, desired_format))
        return super(CachedResource, self).create_response(request, data, response_class, **response_kwargs)

    def wrap_json_response(self, request, data):
        """ Convert a standard JSON response to JSONP
        with a cache while still maintaining the callback param.
        """
        callback = request.GET.get('callback', 'callback')
        return "%s(%s)" % (callback, data)

    def dispatch(self, request_type, request, **kwargs):
        """
        Override of the default dispatch method for a resource.  This
        provides front end level caching for TastyPie APIs.  This requires
        fairly high level tinkering with some API basics but provides
        very effective caching.
        """
        # Since we are hitting this before other validation
        # Ensure that a format has been passed in.

        if getattr(settings, 'API_CACHE_ENABLE', False):

            response_format = request.GET.get('format', None)
            if response_format is None:
                response_format = "json"
            data = self.get_cache(request, force_jsonp=self.is_jsonp(request))
            if data:
                self.is_authenticated(request)
                self.throttle_check(request)
                self.method_check(request, allowed=self._meta.allowed_methods)

                if self.is_jsonp(request):
                    data = self.wrap_json_response(request, data)
                # If the request is a valid cache format, return it.
                return_string = data
                content_type = self._get_valid_content_type(response_format)
                print "cached"
                return HttpResponse(return_string, content_type=content_type, status=200)
            else:
                # If no cache is found go 'Move Along'...
                return super(CachedResource, self).dispatch(request_type, request, **kwargs)
        else:
            return super(CachedResource, self).dispatch(request_type, request, **kwargs)
