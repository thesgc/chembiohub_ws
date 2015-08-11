 cd ~
 
  rm anaconda_requirements.txt pip_requirements.txt

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


