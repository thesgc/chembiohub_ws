from django.conf import settings
from django.contrib.auth.models import User
from django.contrib.sessions.models import Session
from django.http import HttpRequest
from django.utils.importlib import import_module
import datetime


def init_session(session_key):
    """
    Initialize same session as done for ``SessionMiddleware``.
    """
    engine = import_module(settings.SESSION_ENGINE)
    return engine.SessionStore(session_key)


def remove_session_cached_projectlists():
    """
    Remove the session cached project lists from all sessions
    """
    now = datetime.datetime.now()
    request = HttpRequest()

    sessions = Session.objects.filter(expire_date__gt=now)

    for session in sessions:
        username = session.get_decoded().get('_auth_user_id')
        request.session = init_session(session.session_key)
        if "projects_list_cache" in request.session:
            print "found"
            del request.session["projects_list_cache"]
            request.session.save()
