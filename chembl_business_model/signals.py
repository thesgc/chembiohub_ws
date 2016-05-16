__author__ = 'mnowotka'

import django.dispatch

structureChanged = django.dispatch.Signal(providing_args=["instance"])


