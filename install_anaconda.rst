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
  
Now install openbabel and add it to the python path

  wget http://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/openbabel-2.3.2.tar.gz
  
  tar -xvf openbabel-2.3.2.tar.gz
  
  cd openbabel-2.3.2
  
  mkdir build
  
  cd build
  
  cmake .. -DPYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=~/Tools/openbabel-install
  
  make
  
  make install
  
  echo $PYTHONPATH=g >> ~/.bashrc
