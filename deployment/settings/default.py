from .base import *
DEBUG=True
TEMPLATE_DEBUG = DEBUG

import os
import pwd

def get_username():
    return pwd.getpwuid( os.getuid() )[ 0 ]


DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': 'chembiohub_db', # Or path to database file if using sqlite3.
        'USER': get_username(), # Not used with sqlite3.
        'PASSWORD': '', # Not used witis oracle
        'HOST': '/home/chembiohub/postgressocket/', # Set to empty string for localhost. Not used with sqlite3.
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



SESSION_ENGINE = "django.contrib.sessions.backends.cache"
SESSION_CACHE_ALIAS = "default"

RQ_QUEUES = {
    'default': {
        'USE_REDIS_CACHE': 'default',
    },
}


SESSION_COOKIE_NAME = 'chembiohub_sessionid'
CSRF_COOKIE_NAME = 'chembiohubcsrftoken'
STATIC_ROOT = '/srv/chembiohub/chembiohub_ws/deployment/static'
MEDIA_ROOT = '/home/chembiohub/media/'
FLOWJS_PATH = MEDIA_ROOT + 'flow'


LOGIN_REDIRECT_URL = '/chembiohub/#/projects/list'
LOGIN_URL = '/chembiohub/login'
LOGOUT_REDIRECT_URL = '/chembiohub/login'
WEBSERVICES_NAME='chembiohub/api'

INCHI_BINARIES_LOCATION = {"1.02" :"/srv/chembiohub/INCHI-1-BIN/linux/64bit/inchi-1"}

SESSION_CACHE_ALIAS= 'chembiohub'
CACHES = {
"default": {
            "BACKEND": "django_redis.cache.RedisCache",
                                "LOCATION": "redis://127.0.0.1:6379/1",
                                                            "OPTIONS": {
                                                                                                    "CLIENT_CLASS": "django_redis.client.DefaultClient",
                                                                                                                                                    }
                                                                                                                                                },




"chembiohub": {
            "BACKEND": "django_redis.cache.RedisCache",
                    "LOCATION": "redis://127.0.0.1:6379/1",
                            "OPTIONS": {
                                        "CLIENT_CLASS": "django_redis.client.DefaultClient",
                                                }
        }
}

ES_PREFIX = 'chembiohub'


 
STATIC_URL = '/chembiohub/static/'

STATICFILES_DIRS = (
'/srv/chembiohub/chembiohub_ws/src/ng-chem/dist',
)
#Add a template dir so that the html content of index.html can be brought in as a static template when in production so login is handled by Django - see base.py in cbh_chembl_ws_extension (Index() view)
TEMPLATE_DIRS = (
'/srv/chembiohub/chembiohub_ws/src/ng-chem',
)

#Set to 'DEBUG' to view all SQL
DEBUG_SQL = 'INFO'

LOGGING = {
    'version': 1,
    'formatters': {
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s %(process)d %(thread)d %(message)s'
        },
        'simple': {
            'format': '%(levelname)s %(message)s'
        },
    },
    'handlers': {
        # 'sentry': {
        #      'level': 'ERROR',
        #      'class': 'raven.contrib.django.raven_compat.handlers.SentryHandler',
        #  },
        'console': {
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple'
            },
        },
    
    'loggers': {
        'django': {
            'handlers': ['console'],
            'level': DEBUG_SQL,
            'propagate': True,
            },
        },
        
    }


# Set your DSN value
RAVEN_CONFIG = {
     'dsn': 'http://799d9560a5a24a6abc5383e8a4435111:ebc6d747d1654709b812974757213e85@163.1.63.22/2',
 }

for app in INSTALLED_APPS:
    LOGGING["loggers"][app] = {
            'handlers': ['console'],
            'level': 'DEBUG',
            'propagate': True,
            }


