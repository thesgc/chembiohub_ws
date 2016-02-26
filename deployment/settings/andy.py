from .default import *


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
   'deployment',
   'chembl_core_db',
   'chembl_core_model',
  'chembl_business_model', 
         'cbh_core_api',

       'flowjs',
       'cbh_core_model',
       'cbh_datastore_model',
  'cbh_chembl_model_extension',    
  'cbh_chem_api',
  'cbh_chembl_id_generator',
  'solo',
    "django_rq",
            'tastypie',


   )

ROOT_URLCONF = 'deployment.urls_v2'