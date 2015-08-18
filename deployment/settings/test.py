from .development import *


DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': 'tester_cbh_chembl', # Or path to database file if using sqlite3.
        'USER': 'chembl', # Not used with sqlite3.
        'PASSWORD': 'chembl', # Not used witis oracle
        'HOST': '127.0.0.1', # Set to empty string for localhost. Not used with sqlite3.
        'PORT': '5432', # Set to empty string for default. Not used with sqlite3.
        'TEST_NAME' : 'tester_cbh_chembl'
    },
}

ES_PREFIX = "test"