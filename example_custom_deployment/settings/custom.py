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
