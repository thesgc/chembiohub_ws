 cd ~
 
  rm anaconda_requirements.txt pip_requirements.txt
wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda-2.3.0-Linux-x86_64.sh

wget https://raw.githubusercontent.com/thesgc/chembiohub_ws/master/anaconda_requirements.txt
wget https://raw.githubusercontent.com/thesgc/chembiohub_ws/master/pip_requirements.txt

  wget https://raw.githubusercontent.com/thesgc/chembiohub_ws/master/scripts/as_sudoable_user.sh
sh as_sudoable_user.sh
###Then change to that directory and add channels

  cd anaconda/bin
  
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

  ./conda create --yes -m -n chembiohub_ws --file=../../anaconda_requirements.txt
source activate chembiohub_ws
pip install -r ../../pip_requirements.txt

wget https://raw.githubusercontent.com/thesgc/mychembl/master/install_core_libs_Ubuntu.sh

sh install_core_libs_Ubuntu.sh

# operating system, currently only Ubuntu and CentOS are supported:
python -mplatform | grep Ubuntu && export AUX_OS_NAME="Ubuntu" || export AUX_OS_NAME="CentOS"

# rdkit tookit (http://www.rdkit.org/) repository location:
export RDKIT_REPO="https://github.com/rdkit/rdkit"

# rdkit release tag:
export RDKIT_RELEASE="Release_2015_03_1"

# indigo toolkit location
export INDIGO_FILENAME="indigo-python-1.1.11-linux"
export INDIGO_LOCATION="https://dl.dropboxusercontent.com/u/10967207/${INDIGO_FILENAME}.zip"
sudo -H pip install numpy

wget https://raw.githubusercontent.com/thesgc/mychembl/master/rdkit_install.sh

sh rdkit_install.sh

sudo su postgres -c "make installcheck"


  cd ~
  
  wget http://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/openbabel-2.3.2.tar.gz
  
  tar -xvf openbabel-2.3.2.tar.gz
  
  cd openbabel-2.3.2
  
  mkdir build
  
  cd build
  
  cmake .. -DPYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=~/Tools/openbabel-install
  
  #compile with 8 threads for speed
  
  make -j8
  
  make -j8 install
  

  
  
