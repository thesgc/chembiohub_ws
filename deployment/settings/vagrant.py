from .default import *

DEBUG = True
TEMPLATE_DEBUG = DEBUG

ID_PREFIX = "DEV"

WEBSERVICES_NAME='dev'

LOGIN_REDIRECT_URL = "http://localhost:9000/#/projects/list"
LOGOUT_REDIRECT_URL = "login"