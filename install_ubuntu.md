#!/bin/bash

:'Installing ChemBio Hub with Conda

Step 1:

Create a virtual machine - Ubuntu 14.04 64 bit server and a folder for the app to be installed in'


:'Step 2:

Install elasticsearch'

    sudo apt-get install software-properties-common
    wget https://gist.githubusercontent.com/ricardo-rossi/8265589463915837429d/raw/c4a4d1b1cbbbeeecdf90cb496914e50ce2120e2b/ElasticSearch.sh
    sudo sh ElasticSearch.sh 1.7
    sudo update-rc.d elasticsearch defaults

:'Step 3:

Install the other dependencies'

    sudo add-apt-repository ppa:chris-lea/redis-server -y
    sudo apt-get update
    sudo apt-get install -y redis-server  ruby gem ruby-dev unzip fabric git apache2 libapache2-mod-proxy-html  libxml2-dev supervisor nodejs npm tcl8.5 software-properties-common python-software-properties nodejs-legacy
    sudo update-rc.d redis-server defaults

:'Step 4:

Configure the packages'

    sudo a2enmod proxy proxy_http headers expires rewrite
    sudo gem install compass
    sudo npm install -g bower grunt-cli coffee-script

:'Step 5

Run the install script in your home directory where the first argument is the web folder to install in'

    if [ $USER -ne "vagrant" ]
    then
       bash install_linux64.sh 0 Ubuntu
    elif [ $RD_BASE = "/home/chembl/rdkit" ]
    then
       sudo mkdir /srv/chembiohub
       sudo chown -R vagrant /srv/chembiohub
       cd /srv/chembiohub
       git clone  --recursive  git@github.com:thesgc/chembiohub_ws.git
       cd chembiohub_ws
       bash install_linux64.sh chembiohub Ubuntu
    else
       bash install_linux64.sh dev Ubuntu
    fi
