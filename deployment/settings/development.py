from .base import *


SESSION_EXPIRE_AT_BROWSER_CLOSE = False

DEBUG = True
TEMPLATE_DEBUG = DEBUG

ID_PREFIX = "DEV"

WEBSERVICES_NAME='dev'

LOGIN_REDIRECT_URL = "http://localhost:9000/#/projects/list"
LOGOUT_REDIRECT_URL = "login"

# SESSION_REDIS_HOST = 'localhost'
# SESSION_REDIS_PORT = 6379
# SESSION_REDIS_DB = 0
# SESSION_REDIS_PREFIX = 'session'
# SESSION_ENGINE = 'redis_sessions.session'

#
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

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': 'cbh_chembl', # Or path to database file if using sqlite3.
        'USER': 'chembl', # Not used with sqlite3.
        'PASSWORD': 'chembl', # Not used witis oracle
        'HOST': '127.0.0.1', # Set to empty string for localhost. Not used with sqlite3.
        'PORT': '5432', # Set to empty string for default. Not used with sqlite3.
        'TEST_NAME' : 'tester_cbh_chembl'
    },
}
CSRF_COOKIE_NAME = "devcsrftoken"
SESSION_COOKIE_NAME = "dev_sessionid"

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/home/media/media.lawrence.com/static/"
STATIC_ROOT = '/home/vagrant/chembiohub_ws/deployment/static'

# URL prefix for static files.
# Example: "http://media.lawrence.com/static/"
STATICFILES_DIRS = (
# Put strings here, like "/home/html/static" or "C:/www/django/static".
# Always use forward slashes, even on Windows.
# Don't forget to use absolute paths, not relative paths.
'/home/vagrant/chembiohub_ws/src/ng-chem/dist',
)


TEMPLATE_DIRS = [
'/home/vagrant/chembiohub_ws/src/ng-chem/dist',

]
    
TASTYPIE_DOC_API = 'cbh_datastore_ws.urls.api'

# URL prefix for static files.
# Example: "http://media.lawrence.com/static/"
STATIC_URL = '/devapi/r/'

# List of finder classes that know how to find static files in
# various locations.

INCHI_BINARIES_LOCATION = {"1.02" :"/home/vagrant/INCHI-1-BIN/linux/64bit/inchi-1"}

OPEN_BABEL_EXECUTABLE = "/home/vagrant/openbabel-2.3.2/build/bin/babel"


ANONYMOUS_USER_ID = -1

ALLOWED_HOSTS =["testserver"]


INSTALLED_APPS = (    'tastypie_spore_docs',        'devserver',
) + INSTALLED_APPS 


DEVSERVER_MODULES = (
    'devserver.modules.sql.SQLRealTimeModule',
    'devserver.modules.sql.SQLSummaryModule',
    'devserver.modules.profile.ProfileSummaryModule',

    # 'devserver.modules.profile.LineProfilerModule',
)
DEVSERVER_TRUNCATE_SQL = False

INTERNAL_IPS =('0.0.0.0',)
DEVSERVER_DEFAULT_ADDR = '0.0.0.0'
