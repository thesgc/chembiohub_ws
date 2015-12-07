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

source activate $1

pip install -r pip_requirements.txt

git submodule init 
git submodule update
git submodule foreach git checkout master 
git submodule foreach git pull

RANDOM_PORT=$(python generate_port.py)



SUPERVISOR=$"[program:$1_uwsgi]
command=$CONDA_ENV_PATH/bin/uwsgi  --http  :$RANDOM_PORT  --chmod-socket=664  --module=deployment.wsgi
directory=$(pwd)
environment=PATH=$PATH,CONDA_ENV_PATH=$CONDA_ENV_PATH
user=$USER
autorestart=true
redirect_stderr=true" 
printf "$SUPERVISOR" > /tmp/uwsgi
sudo mv /tmp/uwsgi /etc/supervisor/conf.d/$1_uwsgi_supervisor.conf

mkdir $CONDA_ENV_PATH/var/postgressocket
POSTGRES=$"[program:$1_postgresql]
command=$CONDA_ENV_PATH/bin/postgres  -D  $CONDA_ENV_PATH/var/postgresdata  -c  listen_addresses=''  -c  unix_socket_directories=$CONDA_ENV_PATH/var/postgressocket
user=$USER
autorestart=true" 
printf "$POSTGRES" > /tmp/postgres

sudo mv /tmp/postgres /etc/supervisor/conf.d/$1_postgres_supervisor.conf


REDO APACHE
# APACHE="<Directory $(pwd)/deployment/static/>
#  Options Indexes FollowSymLinks
#  AllowOverride None
#  Require all granted
# </Directory>


# <VirtualHost *:80>

# RewriteEngine On
# RewriteCond %{REQUEST_FILENAME} !-f
# RewriteRule ^/$1$ $1/ [L,R=301]
# RewriteRule ^/\$ $1/ [L,R=301]
# RewriteRule ^\$ $1/ [L,R=301]
# ProxyTimeout 300
# ProxyPassMatch ^/$1/((?!#|\s*\$|index\.html|api|admin|login|webauth|webauthlogout).*)\$ !
# AliasMatch ^/$1/static/(.*)\$ $(pwd)/chembiohub_ws/deployment/static/\$1
# AliasMatch ^/$1/((?!#|\s*\$|index\.html).*)$ $(pwd)/chembiohub_ws/deployment/static/\$1
# ProxyPass /$1/ http://127.0.0.1:9090/$1/
# ProxyPassReverse /$1/ http://127.0.0.1:9090/$1/

# </Virtualhost>" 





if [ "$2" == "Ubuntu" ]; then
    sudo service supervisor restart
    printf "$APACHE" > /etc/apache2/sites-available/$1_chembiohub.conf
    sudo a2ensite $1_chembiohub.conf
    sudo service apache2 reload
fi

if [ "$2" == "Centos" ]; then
    sudo service supervisord restart
    printf "$APACHE" > /etc/httpd/conf.d/$1_chembiohub.conf
    sudo /etc/init.d/httpd graceful
fi

initdb -D $CONDA_ENV_PATH/var/postgresdata
createdb -h $CONDA_ENV_PATH/var/postgressocket/ $1_db
psql  -h $CONDA_ENV_PATH/var/postgressocket -c "create extension hstore;create extension rdkit;" $1_db

python manage.py migrate
python manage.py reindex_compounds
python manage.py reindex_datapoint_classifications
python manage.py createsuperuser
python manage.py collectstatic
