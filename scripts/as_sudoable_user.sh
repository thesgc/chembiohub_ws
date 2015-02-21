

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
export POSTGRES_COMMAND="psql template1 -c 'DROP ROLE IF EXISTS $USER; create user $USER with superuser; CREATE EXTENSION IF NOT EXISTS hstore'"
 sudo su postgres -c $POSTGRES_COMMAND
  

###edit pg_hba.conf and add a line for your user 
export ECHO_COMMAND='echo "local all $USER ident" >> /etc/postgresql/9.3/main/pg_hba.conf'
#If version is only 9.1
export ECHO_2='echo "local all $USER ident" >> /etc/postgresql/9.1/main/pg_hba.conf'

  sudo su postgres -c '$ECHO_COMMAND' ||   sudo su postgres -c '$ECHO_2'


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

cd /tmp
wget https://raw.githubusercontent.com/thesgc/chembiohub_ws/master/scripts/as_chembiohub_user.sh

sudo apt-get install -y flex bison build-essential python-numpy cmake python-dev sqlite3 libsqlite3-dev
libboost-dev libboost-python-dev libboost-regex-dev

sudo su chembiohub -c 'bash as_chembiohub_user.sh'

#export DROPCOMMAND='psql template1 -c "CREATE EXTENSION rdkit; DROP ROLE IF EXISTS $USER;" '
#echo $DROPCOMMAND
#  sudo su postgres -c '$DROPCOMMAND'




