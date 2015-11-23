from .base import *
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
        'HOST': '127.0.0.1', # Set to empty string for localhost. Not used with sqlite3.
        'PORT': '5432', # Set to empty string for default. Not used with sqlite3.
    },
}