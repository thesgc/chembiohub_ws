
cd ~

###Now Install the RDKit globally in order to make the database work
  export RDKIT_SOURCE_ROOT=$HOME/rdkit

tar -xvf /var/cache/wget/RDKit_2014_09_2.tar.gz
mv rdkit-Release_2014_09_2 rdkit 
export RDBASE=$HOME/rdkit
export LD_LIBRARY_PATH=$RDBASE/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$RDBASE:$PYTHONPATH
cd rdkit
mkdir build
cd build



#if [$1 = 'travis']
#then
#cmake -DPYTHON_LIBRARY=/home/chembiohub/anaconda/lib/python2.7/config/libpython2.7.a -DPYTHON_INCLUDE_DIR=/home/chembiohub/anaconda/include/python2.7 -DBOOST_ROOT=/home/chembiohub/anaconda ..
#make -j4 install
#else
cmake -DRDK_BUILD_INCHI_SUPPORT=ON -DBOOST_ROOT=/usr/include ..
make -j4 install
#fi

cd $RDBASE/Code/PgSQL/rdkit
make

cd ~
  
  #chmod +x /tmp/anaconda.sh
  
  bash /var/cache/wget/Anaconda-2.1.0-Linux-x86_64.sh -b



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







