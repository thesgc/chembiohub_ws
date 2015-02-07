#!/bin/bash
set -e
#In order to install chembiohub web services and all chembl dependencies on anaconda, run the following:
#===============================

###Now install all of the dependency apt gets in the environment

  wget https://raw.githubusercontent.com/chembl/mychembl/master/install_core_libs.sh

  sh install_core_libs.sh

###Now add a user for the install

  exit

  sudo su postgres
  
  psql postgres
  
  create user vagrant with superuser;
  
  \q
  
  exit
  
###edit pg_hba.conf and add a line for your user 

  cat "local all vagrant ident" > /etc/postgresql/9.3/main/pg_hba.conf

  sudo service postgresql restart
  

###Now Install the RDKit globally in order to make the database work

  wget https://github.com/chembl/mychembl/blob/master/rdkit_install.sh
  
  sh rdkit_install.sh


###Bower and node

  sudo apt-get install nodejs
  
  sudo apt-get install npm
  
  sudo npm install -g bower

  sudo apt-get install nodejs-legacy

###Add a user and switch to it for the non root installs

   
   sudo su

   sudo apt-get install apache2

    sudo useradd -G www-data -s /bin/bash -m chembiohub
   
    sudo su chembiohub


  
###Now install openbabel and indigo and add them to python path

  cd ~
  
  wget http://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/openbabel-2.3.2.tar.gz
  
  tar -xvf openbabel-2.3.2.tar.gz
  
  cd openbabel-2.3.2
  
  mkdir build
  
  cd build
  
  cmake .. -DPYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=~/Tools/openbabel-install
  
  #compile with 8 threads for speed
  
  make -j8
  
  make install
  
###Indigo like this:

  wget https://dl.dropboxusercontent.com/u/10967207/indigo-python-1.1.11-linux.zip

  unzip indigo-python-1.1.11-linux.zip

  rm indigo-python-1.1.11-linux.zip


  echo "export PYTHONPATH=:/var/www/chembiohub_ws:/home/chembiohub/indigo-python-1.1.11-linux:/home/chembiohub/Tools/openbabel-install/lib"  >> ~/.bashrc 
  
  echo 'export DJANGO_SETTINGS_MODULE="deployment.settings.staging"'  >> ~/.bashrc 

###Inchi binaries like this:

  cd ~/Tools
  
  wget http://www.iupac.org/fileadmin/user_upload/publications/e-resources/inchi/1.03/INCHI-1-BIN.zip
  
  unzip INCHI-1-BIN.zip
  
  gunzip INCHI-1-BIN/linux/64bit/inchi-1.gz
  
  chmod +x INCHI-1-BIN/linux/64bit/inchi-1




   cd /var/www
   
   git clone git@github.com:thesgc/chembiohub_ws.git --recursive
   
   cd chembiohub_ws
   
   git submodule init
   
   git submodule update

###Install anaconda locally:

  cd ~
  
  wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
  
  chmod +x miniconda.sh
  
  sh miniconda.sh
  
###Then change to that directory and add channels

  cd miniconda/bin
  
  ./conda config --add channels https://conda.binstar.org/auto
  
  ./conda config --add channels https://conda.binstar.org/auto
  
  ./conda config --add channels https://conda.binstar.org/bcbio
  
  ./conda config --add channels https://conda.binstar.org/ric
  
  ./conda config --add channels https://conda.binstar.org/minadyn
  
  ./conda config --add channels https://conda.binstar.org/pkgw
  
  ./conda config --add channels https://conda.binstar.org/jacksongs
  
  ./conda config --add channels https://conda.binstar.org/mutirri
  
  ./conda config --add channels https://conda.binstar.org/zero323 
    
###Now create a virtualenv using the conda requirements file

  ./conda create --yes python=2.7.6 -m -n chembiohub_ws --file=/var/www/chembiohub_ws/anaconda_requirements.txt




  
###Now ensure that the setting in deployment/settings/base.py matches the location of the inchi binary file - for this install it is:

 ## INCHI_BINARIES_LOCATION = {"1.02" :"/home/chembiohub/Tools/INCHI-1-BIN/linux/64bit/inchi-1"}

###Next we need to link all of our pip packages that are currently subrepos, we can do this by running:

   source ~/miniconda/bin/activate chembiohub_ws
   
   pip install django-cors-headers
   
   cd /var/www/chembiohub_ws/src/chembl_core_db
   
   python setup.py develop
   
   cd ../chembl_core_model/
   
   python setup.py develop
   
   cd ../chembl_webservices/
   
   python setup.py develop
   
   cd ../chembl_business_model/
   
   python setup.py develop
   
   cd ../standardiser/
   
   python setup.py develop
   
   cd ../chembl_beaker/
   
   python setup.py develop
   
   cd ../cbh_chembl_model_extension/
   
   python setup.py develop
   
   cd ../cbh_chembl_ws_extension/
   
   python setup.py develop
   
   cd ../chembl_extras/
   
   python setup.py develop

   cd ../django-flow/
   
   python setup.py develop


###Now we need to link in the ng-chem package as a bower dependency for the front end. This is done by first installing nodejs and bower 


  
###Next go to the folder in src and run bower install

  cd /home/vagrant/chembiohub_ws/src/ng-chem
  
  bower install
  
###We now add this folder to STATICFILES_DIRS to allow it to be served
  
###You can now make changes to ng-chem in src and have them reflect in the static files for the app more generally

###Now create a secret settings file and add a database user for the app

   create user cbh_chembl_usr with password 'xxxxxx';

   create database cbh_chembl_db;

   grant all privileges on  cbh_chembl_db to cbh_chembl_usr;

   grant all privileges on  database cbh_chembl_db to cbh_chembl_usr;
   
###Now migrate the database for the application by running the following:

   source ~/miniconda/bin/activate [YOUR_ENV_NAME]

   python manage.py mysyncdb
   
   python manage.py migrate flowjs

   python manage.py migrate cbh_chembl_model_extension
   
###In order for mysyncdb to work you must have the setting in your settings file:

####   CORE_TABLES_MANAGED = True
   
####   APP_SPECIFIC_TABLES_MANAGED = True
   

   


