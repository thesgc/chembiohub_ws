from deployment.settings.default import *

ROOT_URLCONF = 'example_custom_deployment.urls'

LOGGING["loggers"]["example_custom_deployment"] = {
            'handlers': ['console'],
            'level': 'DEBUG',
            'propagate': True,
            }
#Note that there is no standard password login possible
AUTHENTICATION_BACKENDS = ('example_custom_deployment.backends.CustomAuthBackend',)

ADMIN_USERS = ('bmarsden','ddamerell',)

INSTALLED_APPS = (
      'django.contrib.auth',


       "django_hstore",

    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',

        'grappelli',
    'django.contrib.admin',
      'solo',
   'deployment',
   'chembl_core_db',
   'chembl_core_model',
  'chembl_business_model', 
         'cbh_core_api',

       'flowjs',
       'cbh_core_model',
  'cbh_chembl_model_extension',    
  'cbh_chem_api',
  'cbh_chembl_id_generator',

            'tastypie',

        'django_q',
        'cbh_tests',
          'cbh_datastore_model', #just to get rid of some things
        'cbh_utils',

   )
