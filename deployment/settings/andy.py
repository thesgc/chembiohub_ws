from .default import *

SESSION_EXPIRE_AT_BROWSER_CLOSE = False





ROOT_URLCONF = 'deployment.urls_v2'
WEBSERVICES_NAME="dev/api"
LOGIN_REDIRECT_URL = '/#/projects/list'


REDIS_FOR_ID_GENERATOR = ENV_NAME