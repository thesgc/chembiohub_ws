from deployment.settings.test import *

ROOT_URLCONF = 'example_custom_deployment.urls'

#Note that there is no standard password login possible
AUTHENTICATION_BACKENDS = ('example_custom_deployment.backends.CustomAuthBackend',)

ADMIN_USERS = ('bmarsden','ddamerell',)