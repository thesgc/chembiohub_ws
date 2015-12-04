from .default import *

WEBSERVICES_NAME='dev'
LOGIN_REDIRECT_URL = "http://localhost:9000/#/projects/list"
LOGOUT_REDIRECT_URL = "login"


SESSION_COOKIE_NAME = 'dev_sessionid'
CSRF_COOKIE_NAME = 'devcsrftoken'

DEBUG=True
TEMPLATE_DEBUG=True