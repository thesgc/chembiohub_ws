#!/bin/bash

CONDATEST=$(type conda | grep -c conda)

set -e

export ENV_NAME=$1
OLD_PATH="$PATH"
OPERATING_SYSTEM=$2


if [ "$CONDATEST" != "1" ]
then
    wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-2.4.0-Linux-x86_64.sh 
    if [ "$RD_BASE" == "/home/chembl/rdkit" ]
    then
        bash Anaconda2-2.4.0-Linux-x86_64.sh -b -p /srv/chembiohub/anaconda2
    else
        bash Anaconda2-2.4.0-Linux-x86_64.sh -b    
    fi
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

export CONDA_ENV_PATH=$(conda info | grep "envs dir" | cut -c 25-)/$ENV_NAME
export PATH=$CONDA_ENV_PATH/bin:$OLD_PATH
printf ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>conda"
printf "\n"
printf $CONDA_ENV_PATH
printf "\n"
printf $PATH

pip install -r pip_requirements.txt

git submodule init 
git submodule update
git submodule foreach git checkout master 
git submodule foreach git pull

FOLDER=$(pwd)

cd $CONDA_ENV_PATH/var
wget http://www.inchi-trust.org/download/104/INCHI-1-BIN.zip 
unzip INCHI-1-BIN.zip
gunzip INCHI-1-BIN/linux/64bit/inchi-1.gz
chmod +x INCHI-1-BIN/linux/64bit/inchi-1


cd $FOLDER

mkdir -p $CONDA_ENV_PATH/var/postgressocket
initdb -D $CONDA_ENV_PATH/var/postgresdata





unset LD_LIBRARY_PATH
unset PYTHONPATH
unset RDBASE
USERSUB=$USER

if [ "$USER" == "root" ]
    then
    USERSUB="vagrant"
fi

if [ "$USER" != "travis" ]
 then


    RANDOM_PORT=$(python generate_port.py)


    SUPERVISOR="[program:${ENV_NAME}_uwsgi]
command=$CONDA_ENV_PATH/bin/uwsgi  --http  :$RANDOM_PORT  --chmod-socket=664  --module=deployment.wsgi
directory=$(pwd)
environment=PATH=$PATH,CONDA_ENV_PATH=$CONDA_ENV_PATH
user=$USERSUB
autorestart=true
redirect_stderr=true

[program:${ENV_NAME}_qcluster]
command=$CONDA_ENV_PATH/bin/python manage.py qcluster
directory=$(pwd)
environment=PATH=$PATH,CONDA_ENV_PATH=$CONDA_ENV_PATH
user=$USERSUB
autorestart=true
redirect_stderr=true



" 
    printf "$SUPERVISOR" > /tmp/uwsgi

    POSTGRES="[program:${ENV_NAME}_postgresql]
command=$CONDA_ENV_PATH/bin/postgres  -D  $CONDA_ENV_PATH/var/postgresdata  -c  listen_addresses=''  -c  unix_socket_directories=$CONDA_ENV_PATH/var/postgressocket -c fsync=off -c synchronous_commit=off -c full_page_writes=off -c shared_buffers=1024MB  -c work_mem=256MB
user=$USERSUB
autorestart=true" 
    printf "$POSTGRES" > /tmp/postgres


    sudo mv /tmp/uwsgi "/etc/supervisor/conf.d/${ENV_NAME}_uwsgi_supervisor.conf"

    sudo mv /tmp/postgres "/etc/supervisor/conf.d/${ENV_NAME}_postgres_supervisor.conf"

  #REDO APACHE
    EXCLAM='!'

    APACHE="<VirtualHost *:80>
    <Directory $(pwd)/deployment/static/>
     Options Indexes FollowSymLinks
     AllowOverride None
     Require all granted
    </Directory>

    RewriteEngine On
    RewriteRule ^/$ENV_NAME\$ $ENV_NAME/ [L,R=301]
    RewriteRule ^/\$ $ENV_NAME/ [L,R=301]
    RewriteRule ^\$ $ENV_NAME/ [L,R=301]
    ProxyTimeout 300
    ProxyPassMatch ^/$ENV_NAME/((?${EXCLAM}#|\s*\$|index\.html|api|admin|login|webauth|webauthlogout).*)\$ $EXCLAM
    AliasMatch ^/$ENV_NAME/static/(.*)\$ $(pwd)/deployment/static/\$1
    AliasMatch ^/$ENV_NAME/((?${EXCLAM}#|\s*\$|index\.html).*)\$ $(pwd)/deployment/static/\$1
    ProxyPass /$ENV_NAME/ http://127.0.0.1:$RANDOM_PORT/$ENV_NAME/
    ProxyPassReverse /$ENV_NAME/ http://127.0.0.1:$RANDOM_PORT/$ENV_NAME/

    </VirtualHost>"

fi


#if [ $2 == "Ubuntu" ]
 #then
    #sudo service supervisor restart

#fi

if [ "$OPERATING_SYSTEM" == "Centos" ]
 then
    sudo service supervisord restart

fi




  

if [ "$USER" == "travis" ]
    then

    $CONDA_ENV_PATH/bin/postgres  -D  $CONDA_ENV_PATH/var/postgresdata  -c  listen_addresses=''  -c  unix_socket_directories=$CONDA_ENV_PATH/var/postgressocket &


sleep 5
psql  -h $CONDA_ENV_PATH/var/postgressocket -c "create extension if not exists hstore;create extension if not exists  rdkit;" template1

createdb -h $CONDA_ENV_PATH/var/postgressocket/ ${ENV_NAME}_db -T template1

else
if [ "$OPERATING_SYSTEM" == "Ubuntu" ]
 then
    sudo service supervisor restart
    printf "$APACHE" > /tmp/apache
    sudo mv /tmp/apache /etc/apache2/sites-available/${ENV_NAME}_chembiohub.conf
    sudo a2ensite ${ENV_NAME}_chembiohub.conf
    sudo  a2dissite 000-default
    sudo service apache2 reload
fi




if [  "$OPERATING_SYSTEM"  == "Centos" ]
then
    sudo service supervisord restart
    printf "$APACHE" > /etc/httpd/conf.d/$ENV_NAME_chembiohub.conf
    sudo /etc/init.d/httpd graceful
fi





fi


