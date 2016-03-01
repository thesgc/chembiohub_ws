# coding=utf-8

import datetime
import logging
logger = logging.getLogger('tastypie')

from django.views.debug import CLEANSED_SUBSTITUTE
from django.conf import settings

from tastypie.utils.mime import determine_format
from tastypie.exceptions import UnsupportedFormat
from tastypie.serializers import Serializer


SENSITIVE_PARAMS = [
    'password', 'old_password',]
SENSITIVE_HEADERS = ['Authorization', 'HTTP_AUTHORIZATION']




class ResponseLoggingMiddleware(object):

    def process_response(self, request, response):
        if hasattr(request, 'user') and request.user.is_authenticated():
            credentials = request.user.pk
        else:
            credentials = -1
        serializer = Serializer()
        format = determine_format(
            request, serializer, default_format='application/json')
        # clean POST received from client
        if request.method in ['POST', 'PUT']:
            try:
                cleansed = serializer.deserialize(
                    request.body,
                    format=format)
            except ValueError:
                cleansed = request.POST.dict()
            except UnsupportedFormat:
                cleansed = 'UNSUPPORTED_FORMAT'
            except AttributeError:
                cleansed = 'UNSUPPORTED_STRING'
            except Exception as e:
                logger.warning(
                    u"[RESPONSE] %(c)s: %(e)s",
                    {'c': e.__class__, 'e': e})
                cleansed = 'STRING'
            finally:
                if hasattr(cleansed, 'keys'):
                    cleansed.update(
                        dict(
                            (param, CLEANSED_SUBSTITUTE)
                            for param in SENSITIVE_PARAMS
                            if param in cleansed))
        else:
            cleansed = 'STRING'
        # clean response sent to client
        if response.content and settings.FULL_RESPONSE_LOGGING:
            try:
                content = response.content
                    
            except AttributeError:
                content = 'UNSUPPORTED_STRING'
            except (UnsupportedFormat, ValueError):
                content = 'UNSUPPORTED_FORMAT'
        else:
            content = 'HIDDEN'
        # clean headers from client
        if settings.FULL_HEADER_LOGGING:
                cleansed_headers = getattr(request, 'META', {})
                cleansed_headers.update(
                    dict(
                        (param, CLEANSED_SUBSTITUTE)
                        for param in SENSITIVE_HEADERS
                        if param in cleansed_headers))
        else:
            cleansed_headers = ''

        if 'HTTP_X_FORWARDED_FOR' in request.META:
            ip = request.META['HTTP_X_FORWARDED_FOR'].split(",")[0]
        else:
            ip = request.META['REMOTE_ADDR']

        logger.info(
            u"[ACCESS] %(h)s %(l)s user=%(u)s - %(t)s \"%(r)s\" status=%(s)s "
            u"%(b)s - \"%(f)s\" - \"%(a)s\" - data=%(i)s - headers=%(j)s -"
            u"content=%(o)s",
            {
                'h': ip,
                'l': '-',
                'u': credentials,
                't': datetime.datetime.now(),
                'r': '%s %s' % (request.method, request.get_full_path()),
                's': response.status_code,
                'b': request.META.get('CONTENT_LENGTH', ''),
                'f': request.META.get('HTTP_REFERER', ''),
                'a': request.META.get(
                    'HTTP_USER_AGENT', '').decode('utf-8', 'replace'),
                'i': cleansed,
                'o': content,
                'j': cleansed_headers}
        )
        return response