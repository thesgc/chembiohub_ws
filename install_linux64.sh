#!/bin/bash


if ! hash conda 2>/dev/null; then
    wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-2.4.0-Linux-x86_64.sh
    bash Anaconda2-2.4.0-Linux-x86_64.sh -b    
    rm Anaconda2-2.4.0-Linux-x86_64.sh
     echo "export PATH=$HOME/anaconda2/bin/:\$PATH" >> ~/.bashrc
     source ~/.bashrc
     conda config --add channels https://conda.anaconda.org/jeprescottroy
     conda config --add channels https://conda.anaconda.org/rdkit
     conda config --add channels https://conda.anaconda.org/clyde_fare
fi

conda create --file anaconda_requirements.txt -n $1


source bin/activate $1

pip install -r pip_requirements.txt

git submodule init 
git submodule update
git submodule foreach git checkout master 
git submodule foreach git pull

RANDOM_PORT=${python generate_port.py}



SUPERVISOR="[program:$1_uwsgi]
command=$CONDA_ENV_PATH/bin/uwsgi --http :RANDOM_PORT --chmod-socket=664 --module=deployment.wsgi
directory=${PWD}
environment=PATH=$PATH
user=$USER
autorestart=true
redirect_stderr=true" 
sudo echo $SUPERVISOR > /etc/supervisor/conf.d/$1_uwsgi_supervisor.conf


POSTGRES="[program:$1_postgresql]
command=$CONDA_ENV_PATH/bin/postgres -D $CONDA_ENV_PATH/var/postgresdata -c listen_addresses='' -c unix_socket_directories=$CONDA_ENV_PATH/var/postgressocket
user=$USER
autorestart=true" 
sudo echo $POSTGRES > /etc/supervisor/conf.d/$1_uwsgi_supervisor.conf

sudo supervisorctl reread
sudo supervisorctl reload


initdb -D $CONDA_ENV_PATH/var/postgresdata
createdb -h $CONDA_ENV_PATH/var/postgressocket/ $1_db
psql  -h $CONDA_ENV_PATH/var/postgressocket -c "create extension hstore;create extension rdkit;" $1_db

python manage.py migrate
python manage.py reindex_compounds
python manage.py reindex_datapoint_classifications
python manage.py createsuperuser
python manage.py collectstatic