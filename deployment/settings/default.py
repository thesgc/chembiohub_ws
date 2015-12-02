from .base import *
import os
import pwd

DEBUG=True

TEMPLATE_DEBUG = DEBUG

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))


def get_username():
    return pwd.getpwuid( os.getuid() )[ 0 ]

CONDA_ENV_PATH = os.getenv("CONDA_ENV_PATH")

ENV_NAME = os.path.split(CONDA_ENV_PATH)[1]

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': '%s_db' % ENV_NAME, # Or path to database file if using sqlite3.
        'USER': get_username(), # Not used with sqlite3.
        'PASSWORD': '', # Not used witis oracle
        'HOST': os.getenv("CONDA_ENV_PATH") + '/var/postgressocket', # Set to empty string for localhost. Not used with sqlite3.
        'PORT': '', # Set to empty string for default. Not used with sqlite3.
    },
}


CACHES = {
    "default": {
        "BACKEND": "django_redis.cache.RedisCache",
        "LOCATION": "redis://127.0.0.1:6379/1",
        "OPTIONS": {
            "CLIENT_CLASS": "django_redis.client.DefaultClient",
        }
    }
}



LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
        },
    },
    'loggers': {
        'django': {
            'handlers': ['console'],
            'level': os.getenv('DJANGO_LOG_LEVEL', 'INFO'),
        },
    },
}

SESSION_ENGINE = "django.contrib.sessions.backends.cache"
SESSION_CACHE_ALIAS = "default"

RQ_QUEUES = {
    'default': {
        'USE_REDIS_CACHE': 'default',
    },
}




SESSION_COOKIE_NAME = '%s_sessionid' % ENV_NAME
CSRF_COOKIE_NAME = '%s_csrftoken' % ENV_NAME

STATIC_ROOT = '%s/deployment/static' % BASE_DIR
MEDIA_ROOT = '%s/var/media/' % CONDA_ENV_PATH
FLOWJS_PATH = MEDIA_ROOT + 'flow'


LOGIN_REDIRECT_URL = '/%s/#/projects/list' % ENV_NAME
LOGIN_URL = '/%s/login' %  ENV_NAME
LOGOUT_REDIRECT_URL = '/%s/login' %  ENV_NAME
WEBSERVICES_NAME='%s/api' %  ENV_NAME

INCHI_BINARIES_LOCATION = {"1.02" :"%s/var/INCHI-1-BIN/linux/64bit/inchi-1" % CONDA_ENV_PATH}

SESSION_CACHE_ALIAS= ENV_NAME
CACHES = {
"default": {
            "BACKEND": "django_redis.cache.RedisCache",
                                "LOCATION": "redis://127.0.0.1:6379/1",
                                                            "OPTIONS": {
                                                                                                    "CLIENT_CLASS": "django_redis.client.DefaultClient",
                                                                                                                                                    }
                                                                                                                                                },




ENV_NAME: {
            "BACKEND": "django_redis.cache.RedisCache",
                    "LOCATION": "redis://127.0.0.1:6379/1",
                            "OPTIONS": {
                                        "CLIENT_CLASS": "django_redis.client.DefaultClient",
                                                }
        }
}

ES_PREFIX = ENV_NAME


 
STATIC_URL = '/%s/static/' % ENV_NAME

STATICFILES_DIRS = (
'%s/src/ng-chem/dist' % BASE_DIR,
)
#Add a template dir so that the html content of index.html can be brought in as a static template when in production so login is handled by Django - see base.py in cbh_chembl_ws_extension (Index() view)
TEMPLATE_DIRS = (
'%s/src/ng-chem/' % BASE_DIR,
)

try:
    from .secret import *
except ImportError:
    print "No Secret settings, using default secret key which is insecure"
