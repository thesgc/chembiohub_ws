#!/bin/sh
export PYTHONPATH="/var/www/chembiohub_ws:/home/chembiohub/indigo-python-1.1.11-linux:/home/chembiohub/Tools/openbabel-install/lib" 
exec /home/chembiohub/anaconda/envs/chembiohub_ws/bin/gunicorn  deployment.wsgi:application --debug --log-level debug --preload -b 0.0.0.0:8083  --workers 8 --error-logfile -
