

#!/bin/bash
#set -e
#In order to install chembiohub web services and all chembl dependencies on anaconda, run the following:
#===============================
###First create the user that will run all the code



cd ~
###Now install all of the dependency apt gets in the environment

  wget https://raw.githubusercontent.com/chembl/mychembl/master/install_core_libs.sh

  sed "s/gem install gist//g" install_core_libs.sh >> install_core_libs1.sh
 sh install_core_libs1.sh

###Now add a user for the install
POSTGRES_COMMAND = "psql template1 -c 'DROP ROLE IF EXISTS $USER; create user $USER with superuser; CREATE EXTENSION IF NOT EXISTS hstore'"
 sudo su postgres -c $POSTGRES_COMMAND
  
  
###edit pg_hba.conf and add a line for your user 
ECHO_COMMAND = 'echo "local all $USER ident" >> /etc/postgresql/9.3/main/pg_hba.conf'
  sudo su postgres -c $ECHO_COMMAND

  sudo service postgresql restart
  
   sudo apt-get install -y apache2 
sudo apt-get install -y libapache2-mod-proxy-html libxml2-devsudo apt-get install libxml2-dev
sudo a2enmod proxy proxy_http headers
  sudo apt-get install -y supervisor

  sudo apt-get install -y nodejs

  sudo apt-get install -y npm

  sudo npm install -g bower

  sudo apt-get install -y nodejs-legacy

cd ~
wget http://bitbucket.org/eigen/eigen/get/2.0.15.tar.bz2
tar xvf 2.0.15.tar.bz2
cd eigen-eigen-0938af7840b0/; mkdir build; cd build; cmake ..
sudo make install  # Does the job.
cd ../..
sudo apt-get install pkg-config

id -u chembiohub &>/dev/null ||   sudo useradd -G www-data -s /bin/bash -m chembiohub
   
sudo su chembiohub -c chembiohub.sh

DROP_COMMAND = "psql template1 -c 'DROP ROLE IF EXISTS $USER; '"

  sudo su postgres -c $DROP_COMMAND




