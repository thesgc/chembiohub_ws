#Copyright 2012 Do@. All rights reserved.
#
#Redistribution and use in source and binary forms, with or without modification, are
#permitted provided that the following conditions are met:
#
#   1. Redistributions of source code must retain the above copyright notice, this list of
#      conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above copyright notice, this list
#      of conditions and the following disclaimer in the documentation and/or other materials
#      provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY Do@ ``AS IS'' AND ANY EXPRESS OR IMPLIED
#WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
#FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
#CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
#ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
#ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#The views and conclusions contained in the software and documentation are those of the
#authors and should not be interpreted as representing official policies, either expressed
#or implied, of Do@.

"""
Various Util functinos
"""
__author__ = 'dvirsky'

import redis

'''
Generic singleton wrapper
Created on Sep 21, 2011

@author: dvirsky
'''

import redis
from django.conf import settings 
import django_redis

class Rediston(object):
    """
    This is a base class that holds a connection pool and redis connection instances
    """

    __connPool = None

    _host = 'localhost'
    _port = 6379
    _db = 0
    _timeout = None

    @classmethod
    def _getConnection(cls, mode = 'master'):

        if not hasattr(cls, 'redis'):
            cls.redis = None

        if not cls.redis:
            cls.redis = django_redis.get_redis_connection(settings.REDIS_FOR_ID_GENERATOR)

        return cls.redis

    @classmethod
    def _getPipeline(cls,  mode = 'master', transaction = False):
        """
        Create a pipeline object
        @TODO: support master/slave
        """
        conn = cls._getConnection(mode=mode)
        return conn.pipeline(transaction=transaction)



    @classmethod
    def config(cls, host, port, db, timeout = None):

        cls._host = host
        cls._port = port
        cls._db = db
        cls._timeout = timeout


        cls.__connPool = None


    def resetPool(self):
        """
        hard reconnect to avoid forked sockets etc
        """
        self.__class__.__connPool = None
        self.redis = None
        self.__connect()


    def flush(self):

        if getattr(self, 'pipeline', None):
            self.pipeline.execute()





def InstanceCache(method):
    """
    Use this to wrap instance methods that return a value that doesn't change after the first calculation
    """

    def wrapped(*args, **kwargs):

        _hash = hash((args[1:], tuple(sorted(kwargs.items()))))
        _self = args[0]
        k = '__MethodCache__%s_%x' % (method.__name__, _hash)

        if hasattr(_self, k):
            return getattr(_self, k)

        ret = method(*args, **kwargs)
        setattr(_self, k, ret)
        return ret

    return wrapped

