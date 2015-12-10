#!/bin/bash
ENV_NAME=$1
OLD_PATH="$PATH"
if ! hash conda 2>/dev/null; then
    
    wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-2.4.0-Linux-x86_64.sh 
    bash Anaconda2-2.4.0-Linux-x86_64.sh -b    
    rm Anaconda2-2.4.0-Linux-x86_64.sh
     echo "export PATH=$HOME/anaconda2/bin/:\$PATH" >> ~/.bashrc
     #source ~/.bashrc

     export PATH=$HOME/anaconda2/bin/:$PATH
     conda config --add channels https://conda.anaconda.org/jeprescottroy
     conda config --add channels https://conda.anaconda.org/rdkit
     conda config --add channels https://conda.anaconda.org/clyde_fare
fi

conda create -q -y --file anaconda_requirements.txt -n $ENV_NAME

#source activate $ENV_NAME

export CONDA_ENV_PATH=$HOME/anaconda2/envs/$ENV_NAME
export PATH=$CONDA_ENV_PATH/bin:$OLD_PATH

pip install -r pip_requirements.txt

git submodule init 
git submodule update
git submodule foreach git checkout master 
git submodule foreach git pull




RANDOM_PORT=$(python generate_port.py)


SUPERVISOR="[program:${ENV_NAME}_uwsgi]
command=$CONDA_ENV_PATH/bin/uwsgi  --http  :$RANDOM_PORT  --chmod-socket=664  --module=deployment.wsgi
directory=$(pwd)
environment=PATH=$PATH,CONDA_ENV_PATH=$CONDA_ENV_PATH
user=$USER
autorestart=true
redirect_stderr=true" 
printf "$SUPERVISOR" > /tmp/uwsgi

mkdir $CONDA_ENV_PATH/var/postgressocket
POSTGRES="[program:${ENV_NAME}_postgresql]
command=$CONDA_ENV_PATH/bin/postgres  -D  $CONDA_ENV_PATH/var/postgresdata  -c  listen_addresses=''  -c  unix_socket_directories=$CONDA_ENV_PATH/var/postgressocket
user=$USER
autorestart=true" 
printf "$POSTGRES" > /tmp/postgres


sudo mv /tmp/uwsgi /etc/supervisor/conf.d/$ENV_NAME_uwsgi_supervisor.conf

sudo mv /tmp/postgres /etc/supervisor/conf.d/$ENV_NAME_postgres_supervisor.conf



#REDO APACHE
EXCLAM='!'
APACHE="<Directory $(pwd)/deployment/static/>
 Options Indexes FollowSymLinks
 AllowOverride None
 Require all granted
</Directory>



RewriteEngine On
RewriteCond %{REQUEST_FILENAME} !-f
RewriteRule ^/$ENV_NAME\$ $ENV_NAME/ [L,R=301]
RewriteRule ^/\$ $ENV_NAME/ [L,R=301]
RewriteRule ^\$ $ENV_NAME/ [L,R=301]
ProxyTimeout 300
ProxyPassMatch ^/$ENV_NAME/((?${EXCLAM}#|\s*\$|index\.html|api|admin|login|webauth|webauthlogout).*)\$ $EXCLAM
AliasMatch ^/$ENV_NAME/static/(.*)\$ $(pwd)/chembiohub_ws/deployment/static/$ENV_NAME
AliasMatch ^/$ENV_NAME/((?${EXCLAM}#|\s*\$|index\.html).*)\$ $(pwd)/chembiohub_ws/deployment/static/$ENV_NAME
ProxyPass /$ENV_NAME/ http://127.0.0.1:$RANDOM_PORT/$ENV_NAME/
ProxyPassReverse /$ENV_NAME/ http://127.0.0.1:$RANDOM_PORT/$ENV_NAME/"





if [ "$2" == "Ubuntu" ]; then
    sudo service supervisor restart
    printf "$APACHE" > /etc/apache2/sites-available/$ENV_NAME_chembiohub.conf
    sudo a2ensite $ENV_NAME_chembiohub.conf
    sudo service apache2 reload
fi

if [ "$2" == "Centos" ]; then
    sudo service supervisord restart
    printf "$APACHE" > /etc/httpd/conf.d/$ENV_NAME_chembiohub.conf
    sudo /etc/init.d/httpd graceful
fi

initdb -D $CONDA_ENV_PATH/var/postgresdata
createdb -h $CONDA_ENV_PATH/var/postgressocket/ $ENV_NAME_db
psql  -h $CONDA_ENV_PATH/var/postgressocket -c "create extension hstore;create extension rdkit;" $ENV_NAME_db

python manage.py migrate
python manage.py reindex_compounds
python manage.py reindex_datapoint_classifications
python manage.py createsuperuser
python manage.py collectstatic
