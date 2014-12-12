In order to install chembiohub web services and all chembl dependencies on anaconda, run the following:
===============================

Clone the repository

Install anaconda locally:

  wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
  
  chmod +x miniconda.sh
  
  sh miniconda.sh
  
Then change to that directory and add channels

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
    
Now create a virtualenv using the conda requirements file

  ./conda create --yes python=2.7.6 -m -n beaker --file=/var/www/chembl_beaker/anaconda_requirements.txt

Now install all of the dependency apt gets in the environment

  wget https://raw.githubusercontent.com/chembl/mychembl/master/install_core_libs.sh

  sh install_core_libs.sh

Now add a user for the install

  sudo su postgres
  
  psql postgres
  
  create user astretton with superuser;
  
  \\q
  
  exit
  
edit pg_hba.conf and add a line for your user 

  sudo vim /etc/postgresql/9.3/main/pg_hba.conf
  local all astretton ident

Now Install the RDKit globally in order to make the database work

  wget https://github.com/chembl/mychembl/blob/master/rdkit_install.sh
  
  sh rdkit_install.sh
  
Now install openbabel and indigo and add them to python path

  wget http://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/openbabel-2.3.2.tar.gz
  
  tar -xvf openbabel-2.3.2.tar.gz
  
  cd openbabel-2.3.2
  
  mkdir build
  
  cd build
  
  cmake .. -DPYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=~/Tools/openbabel-install
  
  make
  
  make install
  
Indigo like this:

  wget https://dl.dropboxusercontent.com/u/10967207/indigo-python-1.1.11-linux.zip

  unzip indigo-python-1.1.11-linux.zip

  rm indigo-python-1.1.11-linux.zip

  echo "export PYTHONPATH=:/home/vagrant/chembiohub_ws/indigo-python-1.1.11-linux:/home/vagrant/Tools/openbabel-install/lib" >> ~/.bashrc 

Now we need to link in the ng-chem package as a bower dependency for the front end. This is done by first installing nodejs and bower 


  sudo apt-get install nodejs
  
  sudo apt-get install npm
  
  sudo npm install -g bower

  sudo apt-get install nodejs-legacy
  
Next go to the folder in src and run bower link to create a symlink 

  cd /home/vagrant/chembiohub_ws/src/ng-chem
  
  bower link
  
Next install the symlinked packacke in the static files folder with:

  cd /home/vagrant/chembiohub_ws/deployment/static
  
  bower link ng-chem
  
You can now make changes to ng-chem in src and have them reflect in the static files for the app more generally


