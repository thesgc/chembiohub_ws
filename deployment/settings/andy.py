from .default import *

SESSION_EXPIRE_AT_BROWSER_CLOSE = False

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
            'tastypie',
    'django_q',


   )

Q_CLUSTER = {
    'name': 'DJRedis',
    'workers': 12,
    'timeout': 90,
    'django_redis': ENV_NAME
}

ROOT_URLCONF = 'deployment.urls_v2'