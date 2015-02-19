cd ~
mkdir boost
wget -O boost_1_56_0.tar.gz http://sourceforge.net/projects/boost/files/boost/1.56.0/boost_1_56_0.tar.gz/download
tar xzvf boost_1_56_0.tar.gz
cd boost_1_56_0/
./bootstrap.sh --with-libraries=python,regex --prefix=/home/chembiohub/boost
./bjam install

cd ~

###Now Install the RDKit globally in order to make the database work
  export RDKIT_SOURCE_ROOT=$HOME/rdkit

wget http://sourceforge.net/projects/rdkit/files/rdkit/Q3_2014/RDKit_2014_09_2.tar.gz
tar -xvf RDKit_2014_09_2.tar.gz
mv rdkit-Release_2014_09_2 rdkit 
export RDBASE=$HOME/rdkit
export LD_LIBRARY_PATH=$RDBASE/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$RDBASE:$PYTHONPATH
cd rdkit
mkdir build
cd build
cmake -DRDK_BUILD_INCHI_SUPPORT=ON -DCMAKE_INSTALL_PREFIX=/home/chembiohub/boost/ -DBoost_INCLUDE_DIR=/home/chembiohub/boost/include .. ##-DBOOST_ROOT=/usr/include 
make -j4 installcheck

###Bower and node

  
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
  
  make -j8 install
  
###Indigo like this:

cd ~
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


###Install anaconda locally:

  cd ~
  rm anaconda_requirements.txt pip_requirements.txt

wget https://raw.githubusercontent.com/thesgc/chembiohub_ws/master/anaconda_requirements.txt
wget https://raw.githubusercontent.com/thesgc/chembiohub_ws/master/pip_requirements.txt
  wget http://09c8d0b2229f813c1b93-c95ac804525aac4b6dba79b00b39d1d3.r79.cf1.rackcdn.com/Anaconda-2.1.0-Linux-x86_64.sh -O anaconda.sh
  
  chmod +x anaconda.sh
  
  bash anaconda.sh -b
  
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

  ./conda create --yes python=2.7.6 -m -n chembiohub_ws --file=../../anaconda_requirements.txt
source activate chembiohub_ws
pip install -r ../../pip_requirements.txt




