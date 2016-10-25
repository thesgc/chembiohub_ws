from .default import *

DEBUG = True
TEMPLATE_DEBUG = DEBUG

ID_PREFIX = "DEV"

WEBSERVICES_NAME="dev/api"

ROOT_URLCONF = 'deployment.urls_v2'

LOGIN_REDIRECT_URL = 'http://localhost:9000/#/projects/list'
