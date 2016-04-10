# Django settings for chembl_webservices project.
import os, sys
import threading
threading._DummyThread._Thread__stop = lambda x: 42

TASTYPIE_ALLOW_MISSING_SLASH = True
TASTYPIE_CANNED_ERROR = "An internal server error occurred. Please contact ChEMBL help."


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
  'cbh_chembl_model_extension',    
  'cbh_chem_api',
  'cbh_chembl_id_generator',
  'solo',
            'tastypie',

        'django_q',

   )


ID_PREFIX = "CBH"
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
EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'

EMAIL_HOST_PASSWORD = ''
EMAIL_HOST_USER = ''
EMAIL_PORT = 25
EMAIL_USE_TLS = False
EMAIL_HOST = 'localhost'


# The email address to send on behalf of
SERVER_EMAIL = 'root@localhost'

# If you're u
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
    'django.contrib.messages.context_processors.messages',
    "deployment.context_processors.static_paths",

)
MIDDLEWARE_CLASSES = [
    # 'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    #'cbh_core_ws.middleware.ResponseLoggingMiddleware',
]

MAX_AUTOCOMPLETE_SIZE = 1000

API_CACHE_ENABLE = False
API_CACHE_LENGTH = 900
# 'django.contrib.sessions.middleware.SessionMiddleware',
# 'django.contrib.auth.middleware.AuthenticationMiddleware',

# 'django.middleware.common.CommonMiddleware',
#     'django.middleware.csrf.CsrfViewMiddleware',
# 'django.contrib.messages.middleware.MessageMiddleware',
#     )
FULL_RESPONSE_LOGGING = True
FULL_HEADER_LOGGING = False
ROOT_URLCONF = 'deployment.urls'
INTERNAL_IPS = ('127.0.0.1',)

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'deployment.wsgi.application'

TEMPLATE_DIRS = (
# Put strings here, like "/home/html/django_templates" or "C:/www/django/templates".
# Always use forward slashes, even on Windows.
# Don't forget to use absolute paths, not relative paths.
)



PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))


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

CKEDITOR_CONFIGS = {
    'default': {
        'toolbar': 'full',
        'height': 200,
        'width': 200,
    },
}


CBH_QUERY_TYPES = [
        {
        'name': 'Keyword',
        'value': 'phrase',
        },
        {'name': 'Pick From List',
        'value' : 'pick_from_list'},

        {
            'name': 'Between',
            'value': 'between',
        },
        {
            'name': 'Greater than',
            'value': 'greater_than',
        },
        {
            'name': 'Less than',
            'value': 'less_than',
        },
        {
            'name': 'Blank',
            'value': 'blanks',
        },
        {
            'name': 'Not blank',
            'value': 'nonblanks',
        }

    ]

CBH_SORT_DIRECTIONS = [{'name': 'No Sort', 'value': "No Sort" },
             {'name': 'Asc', 'value': 'asc'},
             {'name': 'Desc', 'value': 'desc'},]

CBH_HIDE_SHOW = [{'name' : 'Show', 'value': 'show'},
                {'name' : 'Hide & Close', 'value': 'hide'},]

CBH_SORT_SCHEMAFORM = {
    "default" :{
        "form" : [
                    {
                          "key" : "sort_direction",
                          "title": "Sort direction",
                          "type": "select",
                          "titleMap": CBH_SORT_DIRECTIONS,
                          "htmlClass": "col-sm-12",
                          "type": "radiobuttons",
                          "style": {
                            "selected": "btn-success btn-sm",
                            "unselected": "btn-default btn-sm"
                          },
                          "onChange": "sortChanged(modelValue,form)",
                          "disableSuccessState":True,
                          "feedback": False,
                    },            

                ],

        "schema": {
            "type": "object",
            "properties": {
                "field_path" :{
                  "type" : "string",
                  "default" : ""
                },
                "sort_direction" : {
                    "type": "string",
                    "default": "No Sort",
                    "enum" : [sd["value"] for sd in CBH_SORT_DIRECTIONS],
                    },
              
            },
            "required" : [ "sort_direction"]
        }
    }
}





CBH_HIDE_SCHEMAFORM = {
    "default" :{
        "form" : [
                    
                    {
                          "key" : "hide",
                          "type": "select",
                          "titleMap": CBH_HIDE_SHOW,
                          "type": "radiobuttons",
                          "htmlClass": "col-sm-12",
                          "onChange": "hideChanged(modelValue,form)",
                          "disableSuccessState":True,
                          "feedback": False,
                          "style": {
                            "selected": "btn-success btn-sm",
                            "unselected": "btn-default btn-sm"
                          },
                    },
                   

                ],

        "schema": {
            "type": "object",
            "properties": {
                
                "hide" : {
                    'type': 'string',
                    "enum" : [h["value"] for h in CBH_HIDE_SHOW],
                    "title": "Hide / show Column",
                    "default": "show"
                },
               
            },
            "required" : [ ]
        }
    }
}
        

CBH_CHEMICAL_QUERY_TYPES = [
  
        {
            'name': 'Substructure',
            'value': 'with_substructure',
        },
        {
            'name': 'Exact Match',
            'value': 'flexmatch',
        },
]


CBH_CHEMICAL_QUERY_SCHEMAFORM = {
  'default': {
    "form": [

                

                    {
                          "type": "section",
                          "htmlClass": "row",
                          "items": [
                              {
                                    "key" : "query_type",
                                    "title": "Filter type",
                                    "type": "select",
                                    "titleMap": CBH_CHEMICAL_QUERY_TYPES,
                                    "type": "radiobuttons",
                                    "style": {
                                      "selected": "btn-success btn-sm",
                                      "unselected": "btn-default btn-sm"
                                    },
                                    "htmlClass": "col-sm-4",
                                    "onChange": "chemicalUpdated(modelValue,form)",
                                    "disableSuccessState":True,
                                    "feedback": False,
                              },
                              { 
                                'title' : 'Sketch',
                                'key' : 'molfile',
                                'condition' : 'model.query_type=="with_substructure" || model.query_type=="flexmatch"', 
                                "htmlClass": "col-sm-8",
                                "onChange": "not avaliable",
                                "disableSuccessState":True,
                                "feedback": False,
                                 "type" : "chemdoodle",
                              }
                          ]
                    },
                    {
                          "type": "section",
                          "htmlClass": "row",
                          "items": [
                            {
                              "type": "template",
                              "template": '<div class="help-block col-xs-12 has-error">Invalid molecule, check number of bonds around each atom.</div>',
                              "condition" : '!model.inprogress && model.error && model.molfile!=""'
                            },
                            {
                              "type": "template",
                              "template": '<div class="help-block col-xs-12 has-warning">Structure search has changed, Click Apply Filter to search or clear to cancel.</div>',
                              "condition" : 'model.molfile && !model.filter_is_applied && !model.error'
                            },
                            {
                              "type": "template",
                              "template": '<div class="help-block col-xs-12 ">Data below is filtered for this molecule.</div>',
                              "condition" : '!model.inprogress && model.molfile!="" && model.filter_is_applied'
                            },
                            {
                              "type": "template",
                              "template": '<div class="help-block col-xs-12 has-error">Draw a molecule and click Apply Filter to run structure search.</div>',
                              "condition" : '!model.inprogress && model.molfile==""'
                            },
                            
                          ]
                      },
                    {
                          "type": "section",
                          "htmlClass": "row",
                          "items": [
                          {
                                  "type": 'button',
                                  "value": "submit", 
                                  "htmlClass": "col-xs-3 pull-right",
                                  "style": 'btn-success btn-block ', 
                                  "title": 'Apply and Close', 
                                  "condition" : 'model.smiles && !model.inprogress',
                                  "onClick" : "chemicalsubmit();closeMenu()"

                            },
                            {
                              "type": "template",
                              "template": '<div class="col-xs-3 pull-right"><button class="btn btn-success  btn-block " type="button" disabled>Apply and Close</button></div>',
                              "condition" : '!model.smiles || model.inprogress'
                            },
                            {
                                  "type": 'button', 
                                  "htmlClass": "col-xs-3 pull-right",
                                  "style": 'btn-primary btn-block ',
                                  "title": 'Close', 
                                  "onClick": "closeMenu()" 

                            },
                            {
                                  "type": 'button', 
                                  "htmlClass": "col-xs-3 pull-right",
                                  "style": 'btn-danger btn-block ', 
                                  "title": 'Clear', 
                                  "onClick": "removeStructureSearch()" 

                            },
                            
                          ]
                      }
    ],
    'schema': {
          "type": "object",
          "properties": {
          "query_type": {
                        "type": "string",
                        "enum" : [qt["value"] for qt in CBH_CHEMICAL_QUERY_TYPES],
                        "default": "with_substructure"
                    },
          "molfile" : {
                        'type': 'string',
                        'default': ""
                    }
          }
    },
    "required": ["molfile"]
  }
}



CBH_QUERY_SCHEMAFORM = {
    "default" :{
        "form" : [
                    
                    {
                          "key" : "query_type",
                          "title": "Filter type",
                          "type": "select",
                          "titleMap": CBH_QUERY_TYPES,
                          "type": "radiobuttons",
                          "style": {
                            "selected": "btn-success btn-sm",
                            "unselected": "btn-default btn-sm"
                          },
                          "htmlClass": "row",
                          "onChange": "queryTypeChanged(modelValue,form)",
                          "disableSuccessState":True,
                          "feedback": False,
                    },
                    
                    {
                          "type": "section",
                          "htmlClass": "row",
                          "items": [
                    
                    {
                        'title': 'Keyword or phrase',
                        'key': 'phrase',
                        'condition': 'model.query_type=="phrase"',
                        "htmlClass": "col-sm-4",
                        "onChange": "updated(modelValue,form)",
                          "disableSuccessState":True,
                          "feedback": False,

                    },
                   
                    {
                        'key': 'greater_than',
                        'condition': 'model.query_type=="greater_than" || model.query_type=="between"',
                          "htmlClass": "col-sm-4",
                          "onChange": "updated(modelValue,form)",
                          "disableSuccessState":True,
                          "feedback": False,
                    },
                    {
                        'key': 'less_than',
                        'condition': 'model.query_type=="less_than" || model.query_type=="between"',
                        "htmlClass": "col-sm-4",
                        "onChange": "updated(modelValue,form)",
                          "disableSuccessState":True,
                          "feedback": False,
                    }, 
                    { 
                        'title' : 'Pick from list',
                        'key' : 'pick_from_list',
                        'condition' : 'model.query_type=="pick_from_list"', 
                        "htmlClass": "col-sm-6",
                        "onChange": "updated(modelValue,form)",
                        "disableSuccessState":True,
                        "feedback": False,
                        "type" : "filtereddropdown",
                        "options":{
                          "tagging" : False,
                          "fetchDataEventName" : "openedSearchDropdown",
                          "dataArrivesEventName" : "autoCompleteData",
                          "multiple" : True,
                          "staticItems" : [] 
                        }
                    },
                    {
                          "type": "section",
                          "htmlClass": "row pull-right",
                          "items": [

                            {
                                  "type": 'button', 
                                  "htmlClass": "col-xs-6",
                                  "style": 'btn-primary btn-block',
                                  "title": 'Close', 
                                  "onClick": "closeMenu()" 

                            },
                            {
                                  "type": 'button', 
                                  "htmlClass": "col-xs-6",
                                  "style": 'btn-danger btn-block', 
                                  "title": 'Clear', 
                                  "onClick": "clearFilters()" 

                            },
                          ]
                      }
                          
                        ]},
                        

                    

                ],

        "schema": {
            "type": "object",
            "properties": {
                "field_path" :{
                  "type" : "string",
                  "default" : ""
                },
                "pick_from_list":{
                  "type" : "array",
                  "items": {
                    "type": "string"
                  },
                  "default" : []

                },
                "sort_direction" : {
                    "type": "string",
                    "default": "No Sort",
                    "enum" : [sd["value"] for sd in CBH_SORT_DIRECTIONS],
                    },
                "query_type": {
                    "type": "string",
                    "enum" : [qt["value"] for qt in CBH_QUERY_TYPES],
                    "default": "phrase"
                },
                "phrase" : {
                        'type': 'string',
                        "default": ""
                },
                "any_of" : {
                        'type': 'array',
                        "items": {
                          "type": "string",
                          "title" : " ",
                           "disableSuccessState":True,
                  "feedback": False,
                  "required": True,
                  "default": ""
                        }
                },
                "equals" : {
                    'type': 'string',
                    'title' : 'Exactly equal to',
                    "default": ""
                },
                "starts_with" : {
                    'type': 'string',
                    'title' : 'Starts with',
                    "default": ""
                },
                "ends_with" : {
                    'type': 'string',
                    'title' : 'Ends with',
                    "default": ""
                },
                "greater_than" : {
                    'type': 'string',
                    'title' : 'Greater than',
                    "default": ""
                },
                "less_than" : {
                    'type': 'string',
                    'title' : 'Less than',
                    "default": ""
                }
                
            },
            "required" : [ "sort_direction", "query_type", "phrase", "any_of", "equals", "starts_with", "ends_with", "greater_than", "less_than"]
        }
    }
}
        






TABULAR_DATA_SETTINGS = { 
    "seach_page_edit_mode": {
        "start": ["project_counter","properties.archived", "clone", "image","uuid", "projectfull.name"],
        "end": []
    },
    "export" : {
        "start" :["project_counter","image","uuid", "projectfull.name"],
        "end" : ["created_by" ,"timestamp" , "id" , "multiple_batch_id"]
                },
    "cbh.restoreitems": {
        "start" :["project_counter","properties.archived", "image","uuid", "projectfull.name"],
        "end" : ["created_by" ,"timestamp" , "id" , "multiple_batch_id"]
                },
    "cbh.archiveitems": {
        "start" :["project_counter","properties.archived", "image","uuid", "projectfull.name"],
        "end" : ["created_by" ,"timestamp" , "id" , "multiple_batch_id"]
                },
    "cbh.searchv2": {
        "start" :[ "project_counter","clone", "image","uuid", "projectfull.name"],
        "end" : ["created_by" ,"timestamp" , "id" , "multiple_batch_id"]
                },
    "add_page" : {
        "start" :["project_counter","image", "id", "originalSmiles","properties.action","standardInchiKey"],
        "end" : []
                },
    "indexing" : {
        "start" :["project_counter","uuid", "projectfull.name", "properties.archived","projectfull.project_type.saved_search_project_type"],
        "end" : ["created_by" ,"timestamp" , "id" , "multiple_batch_id"]
    },
    "indexing_temp" : {
        "start" :[ "row", "originalSmiles","properties.action","standardInchiKey", ],
        "end" : ["created_by" ,"timestamp" , "id" , "multiple_batch_id"]
    },
    "schema": {
        "project_counter":{
            "knownBy" : "ID in project",
            "data" : "project_counter",
            "editable": False,
            "className": "htCenter htMiddle ",
        },
        "projectfull.project_type.saved_search_project_type" : {
            "noSort": True,
            "knownBy": "Is a saved search",
            "data": "projectfull.project_type.saved_search_project_type",
            "searchFormType" : "pick_from_list",
            "renderer_named": "archivedRenderer",
            "editable": False,
            "className": "htCenter htMiddle ",
        },
        "properties.archived" : {
            "noSort": True,
            "knownBy": "Archive",
            "data": "properties.archived",
            "searchFormType" : "pick_from_list",
            "renderer_named": "archivedRenderer",
            "editable": True,
            "className": "htCenter htMiddle ",
        },
        "clone": {
            "noSort": True,
            "knownBy": "Clone",
            "data": "clone",
            "searchFormType" : None,
            "renderer_named": "cloneRenderer",
            "editable": False,
            "className": "htCenter htMiddle ",
        },
        "image" : {
            "noSort": True,
            "knownBy": "Structure",
            "data": "image",
            "searchFormType" : "chemical",
            "renderer_named": "coverRenderer",
            "editable": False,
            "className": "htCenter htMiddle ",
            "field_type": "b64png"
        },
        "id" : {
            "sortOrder": "none",
            "knownBy": "Batch ID",
            "data": "id",
            "field_type"  : "integer",
            "searchFormType" : "pick_from_list",
            "editable": False,
            "className": "htCenter htMiddle " 
        }, 
        "originalSmiles" : {
            "noSort": True,
            "editable": False,
            "knownBy": "Info",
            "searchFormType" : "upload_info",
            "data": "originalSmiles",
            "renderer_named": "infoRenderer"
        }, 
        "properties.action" : {
            "sortOrder": "none",
            "knownBy": "Action",
            "data": "properties.action",
            "searchFormType" : "pick_from_list",
            "type": "dropdown",
            "editable": True,
            "source": ["New Batch", "Ignore"],
            "className": "htCenter htMiddle "
        }, 
        "standardInchiKey" : {
            "sortOrder": "none",
            "knownBy": "Inchi Key",
            "data": "standardInchiKey",
            "searchFormType" : "pick_from_list",
            "editable": False,
            "renderer_named": "linkRenderer"
        },
        "uuid" : {
            "noSort": True,
            "knownBy": ID_PREFIX + "_ID",
            "data": "uuid",
            "searchFormType" : "pick_from_list",
            "renderer_named": "modalLinkRenderer",
            "editable": False,
            "className": " htCenter htMiddle ",
        },
        "projectfull.name": {
            "knownBy": "Project",
            "data": "projectfull.name",
            "searchFormType" : "default",
            "editable": False,
            "className": "htCenter htMiddle ",
            "renderer_named": "projectRenderer",
        },
        "created_by" : {
            "noSort": True,
            "knownBy": "Added By",
            "data": "created_by",
            "searchFormType" : "pick_from_list",
            "editable": False,
            "className": "htCenter htMiddle ",
        }, 
        "timestamp" : {
            "noSort": True,
            "knownBy": "Date",
            "data": "timestamp",
            "field_type"  : "date",
            "searchFormType" : "date_range",
            "editable": False,
            "className": "htCenter htMiddle ",
        }, 
        "multiple_batch_id" :{
            "noSort": True,
            "knownBy": "Upload ID",
            "data": "multiple_batch_id",
            "searchFormType" : "pick_from_list",
            "editable": False,
            "field_type"  : "integer",
            "className": "htCenter htMiddle ",
        }
    }
}

ELASTICSEARCH_MAX_FIELD_LENGTH = 256

ELASTICSEARCH_INDEX_MAPPING = {
        "settings": {
            "index.store.type": "niofs",
            "analysis" : {
                    "char_filter" : {
                        "special_char_space_out" :{ # put spaces around special characters so they can still be indexed
                            "type":"pattern_replace",
                            "pattern":"([()\[\].,\-\+])",
                            "replacement":" $1 "
                        }
                    },
                    "analyzer" : {
                        "default_index" : {
                            "tokenizer" : "whitespace",
                            "filter" : [
                                "lowercase"
                            ],
                            "char_filter" : [
                                "html_strip", "special_char_space_out"
                            ]
                        }
                    }
                },
        },
        "mappings": {
            "newbatches": {
                "dynamic": False,
                "_all": {"enabled": False},
                "date_detection": False,
                "properties":{
                    "indexed_fields" :{
                        "type": "nested",
                        "properties" : {
                                "name": {
                                    "type": "string", 
                                    "index": "not_analyzed"
                                },
                                "field_path": {
                                    "type": "string", 
                                    "index": "not_analyzed"
                                },
                          
                                "value":  
                                      {
                                        "type": "string", 
                                        "store": "no", 
                                        "index_options": "positions", 
                                        "index": "analyzed", 
                                        "omit_norms": True, 
                                        "analyzer" : "default_index",
                                        "fields": {
                                            "raw": {"type": "string", "store": "no", "index": "not_analyzed", "ignore_above": ELASTICSEARCH_MAX_FIELD_LENGTH}
                                        }
                                    }
                            }
                        }
                },
             "dynamic_templates": [{
                    "ignored_fields": {
                        "match": "*",
                        "match_mapping_type": "string",
                        "mapping": {
                            "type": "string", "store": "no", "include_in_all": False, "index" : "no"
                        }
                    }
                }]
                
            }
        }
    }
