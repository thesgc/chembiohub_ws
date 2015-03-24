

#!/bin/bash
#set -e
#In order to install chembiohub web services and all chembl dependencies on anaconda, run the following:
#===============================
###First create the user that will run all the code
sudo apt-get install -y ruby gem ruby-dev unzip fabric
 sudo gem install compass
sudo mkdir -p /var/cache/wget
sudo chmod ugo+rw -R /var/cache/wget
cd ~
###Now install all of the dependency apt gets in the environment

  wget https://raw.githubusercontent.com/chembl/mychembl/master/install_core_libs.sh

  sed "s/gem install gist//g" install_core_libs.sh >> install_core_libs1.sh
 sh install_core_libs1.sh

###Now add a user for the install

   sudo apt-get install -y apache2 
sudo apt-get install -y libapache2-mod-proxy-html libxml2-devsudo apt-get install libxml2-dev
sudo a2enmod proxy proxy_http headers expires
  sudo apt-get install -y supervisor

  sudo apt-get install -y nodejs

  sudo apt-get install -y npm


  sudo apt-get install -y nodejs-legacy
  sudo apt-get install -y ruby gem ruby-dev
  sudo gem install compass

  sudo npm install -g bower grunt-cli coffee-script

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
sudo wget http://sourceforge.net/projects/rdkit/files/rdkit/Q3_2014/RDKit_2014_09_2.tar.gz -N -P /var/cache/wget

sudo  wget http://09c8d0b2229f813c1b93-c95ac804525aac4b6dba79b00b39d1d3.r79.cf1.rackcdn.com/Anaconda-2.1.0-Linux-x86_64.sh -N -P /var/cache/wget 

#sudo apt-get install -f -y flex bison build-essential python-numpy cmake python-dev sqlite3 libsqlite3-dev
export COMM="bash as_chembiohub_user.sh $USER"
sudo su chembiohub -c "$COMM" 

#export DROPCOMMAND='psql template1 -c "CREATE EXTENSION rdkit; DROP ROLE IF EXISTS $USER;" '
#echo $DROPCOMMAND
#  sudo su postgres -c '$DROPCOMMAND'

export RDBASE=/home/chembiohub/rdkit
export LD_LIBRARY_PATH=$RDBASE/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$RDBASE:$PYTHONPATH
cd /home/chembiohub/rdkit/Code/PgSQL/rdkit/
sudo make install

export POSTGRES_COMMAND="psql template1 -c ' CREATE EXTENSION IF NOT EXISTS hstore; CREATE EXTENSION IF NOT EXISTS rdkit;'"
 sudo su postgres -c "$POSTGRES_COMMAND"

sudo mkdir /var/www/automated_reg_installs
sudo chown -R chembiohub:www-data /var/www/automated_reg_installs
