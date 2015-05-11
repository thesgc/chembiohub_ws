# Django settings for chembl_webservices project.
import os, sys
import threading
threading._DummyThread._Thread__stop = lambda x: 42

TASTYPIE_ALLOW_MISSING_SLASH = True
TASTYPIE_CANNED_ERROR = "An internal server error occurred. Please contact ChEMBL help."

ID_PREFIX = "UOX"
BEAKER_PATH =  "http://staging.chembiohub.ox.ac.uk/utils/"
SESSION_COOKIE_HTTPONLY = False
SESSION_SAVE_EVERY_REQUEST = True
OPEN_SOURCE = True
LOGIN_REDIRECT_URL = '/r/#/projects/list'
LOGOUT_REDIRECT_URL = "login"

DEBUG = False
TEMPLATE_DEBUG = DEBUG

WS_BASE_URL='/chemblws'
WS_DOCS_TITLE='Chem Bio Hub web services based upon ChEMBL web services live documentation'

ADMINS = (
    ('Andrew Stretton', 'andrew.stretton@sgc.ox.ac.uk'),
    )

MANAGERS = ADMINS

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': 'cbh_chembl', # Or path to database file if using sqlite3.
        'USER': 'chembl', # Not used with sqlite3.
        'PASSWORD': 'chembl', # Not used witis oracle
        'HOST': '127.0.0.1', # Set to empty string for localhost. Not used with sqlite3.
        'PORT': '5432', # Set to empty string for default. Not used with sqlite3.
    },
}

EXPORT_MODE = True
CORE_TABLES_MANAGED = True
APP_SPECIFIC_TABLES_MANAGED = True
#COMPOUND_MOLS_TABLE = 'mols_rdkit'
COMPOUND_MOLS_TABLE = 'compound_mols'

#CTAB_COLUMN = 'm'
CTAB_COLUMN = 'ctab'

#CHEMBL_SCHEMA_NAME = u'CHEMBL16_SNAPSHOT' #u'CHEMBL_TMP'

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# In a Windows environment this must be set to your system time zone.
TIME_ZONE = 'America/Chicago'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale.
USE_L10N = True

# If you set this to False, Django will not use timezone-aware datetimes.
USE_TZ = True

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/home/media/media.lawrence.com/media/"
MEDIA_ROOT = '/var/data/chembiohub_ws/'
FLOWJS_PATH = MEDIA_ROOT + 'flow'
# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://media.lawrence.com/media/", "http://example.com/media/"
MEDIA_URL = ''

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/home/media/media.lawrence.com/static/"

# URL prefix for static files.
# Example: "http://media.lawrence.com/static/"
STATIC_URL = '/reg/'

# Additional locations of static files


# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    )

# Make this unique, and don't share it with anybody.
SECRET_KEY = '3v2xb&amp;@&amp;_kibf0o!4m249njy3!qjxptht0m%q2w&amp;ry8v&amp;ok$na'

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
    'django.template.loaders.eggs.Loader',
    )

TEMPLATE_CONTEXT_PROCESSORS = (
    'django.core.context_processors.csrf',
    'django.contrib.auth.context_processors.auth',
    'django.core.context_processors.debug',
    'django.core.context_processors.i18n',
    'django.core.context_processors.media',
    'django.core.context_processors.static',
    'django.core.context_processors.request',
    'django.contrib.messages.context_processors.messages'
)

MIDDLEWARE_CLASSES = (
'django.contrib.sessions.middleware.SessionMiddleware',
'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
'django.contrib.auth.middleware.AuthenticationMiddleware',
'django.contrib.messages.middleware.MessageMiddleware',
    )

ROOT_URLCONF = 'deployment.urls'
INTERNAL_IPS = ('127.0.0.1',)

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'deployment.wsgi.application'

TEMPLATE_DIRS = (
# Put strings here, like "/home/html/django_templates" or "C:/www/django/templates".
# Always use forward slashes, even on Windows.
# Don't forget to use absolute paths, not relative paths.
)
INSTALLED_APPS = (

       "django_hstore",

    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
        'grappelli',

    'django.contrib.admin',
    'tastypie',
   'deployment',


   'chembl_core_db',
   'chembl_core_model',
  'chembl_business_model', 
       'flowjs',


  'cbh_chembl_model_extension',    
  'cbh_chembl_ws_extension',

   )



PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))

# LOGGING = {
#     'version': 1,
#     'disable_existing_loggers': True,
#     'root': {
#         'level': 'WARNING',
#         'handlers': ['sentry'],
#     },
#     'formatters': {
#         'verbose': {
#             'format': '%(levelname)s %(asctime)s %(module)s %(process)d %(thread)d %(message)s'
#         },
#     },
#     'handlers': {
#         'sentry': {
#             'level': 'ERROR',
#             'class': 'raven.contrib.django.raven_compat.handlers.SentryHandler',
#         },
#         'console': {
#             'level': 'DEBUG',
#             'class': 'logging.StreamHandler',
#             'formatter': 'verbose'
#         }
#     },
#     'loggers': {
#         'django.db.backends': {
#             'level': 'ERROR',
#             'handlers': ['console'],
#             'propagate': False,
#         },
#         'raven': {
#             'level': 'DEBUG',
#             'handlers': ['console'],
#             'propagate': False,
#         },
#         'sentry.errors': {
#             'level': 'DEBUG',
#             'handlers': ['console'],
#             'propagate': False,
#         },
#     },
# }
# Set your DSN value
RAVEN_CONFIG = {
    'dsn': 'http://799d9560a5a24a6abc5383e8a4435111:ebc6d747d1654709b812974757213e85@163.1.63.22/2',
}

# Add raven to the list of installed apps
INSTALLED_APPS = INSTALLED_APPS + (
    # ...
    'raven.contrib.django.raven_compat',
)

CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
        'LOCATION': 'ws_cache',
        }
}

CACHE_MIDDLEWARE_SECONDS = 3000000
#import cbh_chembl_ws_extension
#import cbh_chembl_model_extension
INCHI_BINARIES_LOCATION = {"1.02" :"/home/chembiohub/Tools/INCHI-1-BIN/linux/64bit/inchi-1"}

OPEN_BABEL_EXECUTABLE = "/home/chembiohub/openbabel-2.3.2/build/bin/babel"

ALLOWED_INCLUDE_ROOTS = ("/home/vagrant/",'/var/www/chembiohub_ws/')

ALLOWED_HOSTS = [ "staging.chembiohub.ox.ac.uk", "chembiohub.ox.ac.uk"]

#email config
EMAIL_HOST = 'smtp.jr2.ox.ac.uk'
EMAIL_HOST_USER = ''
EMAIL_HOST_PASSWORD = ''
EMAIL_PORT = 587
EMAIL_USE_TLS = True
