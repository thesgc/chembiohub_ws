from deployment.settings.default import *

ROOT_URLCONF = 'example_custom_deployment.urls'

LOGGING["loggers"]["example_custom_deployment"] = {
            'handlers': ['console'],
            'level': 'DEBUG',
            'propagate': True,
            }
