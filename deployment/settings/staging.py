from .base import *

DEBUG = True
TEMPLATE_DEBUG = DEBUG

LOGIN_REDIRECT_URL = '/reg/app/index.html'

ID_PREFIX = "UOXSTA"


# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/home/media/media.lawrence.com/static/"
STATIC_ROOT = '/var/www/chembiohub_ws/deployment/static'

# URL prefix for static files.
# Example: "http://media.lawrence.com/static/"
STATIC_URL = '/reg/'

# Additional locations of static files
STATICFILES_DIRS = (
# Put strings here, like "/home/html/static" or "C:/www/django/static".
# Always use forward slashes, even on Windows.
# Don't forget to use absolute paths, not relative paths.
'/var/www/chembiohub_ws/src/ng-chem',
)

# List of finder classes that know how to find static files in
# various locations.

from .stagingsecret import *
