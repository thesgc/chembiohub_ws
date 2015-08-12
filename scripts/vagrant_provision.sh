sudo apt-get install -y ruby gem ruby-dev unzip fabric
 sudo gem install compass
sudo mkdir -p /var/cache/wget
sudo chmod ugo+rw -R /var/cache/wget
cd ~
###Now install all of the dependency apt gets in the environment
#if [[ "vagrant" != $USER ]]
#  then
#  wget https://raw.githubusercontent.com/chembl/mychembl/master/install_core_libs_Ubuntu.sh

#  sed "s/gem install gist//g" install_core_libs_Ubuntu.sh >> install_core_libs1.sh
# sh install_core_libs1.sh

###Now add a user for the install

 #  sudo apt-get install -y apache2 postgresql-client-common
# fi
#sudo apt-get install -y libapache2-mod-proxy-html libxml2-devsudo apt-get install libxml2-dev
#sudo a2enmod proxy proxy_http headers expires rewrite
#  sudo apt-get install -y supervisor

#  sudo apt-get install -y nodejs

#  sudo apt-get install -y npm



cd ~/Tools
  
  wget http://www.iupac.org/fileadmin/user_upload/publications/e-resources/inchi/1.03/INCHI-1-BIN.zip
  
  unzip INCHI-1-BIN.zip
  
  gunzip INCHI-1-BIN/linux/64bit/inchi-1.gz
  
  chmod +x INCHI-1-BIN/linux/64bit/inchi-1

sudo apt-get install openjdk-7-jre-headless -y
wget -qO - https://packages.elastic.co/GPG-KEY-elasticsearch | sudo apt-key add -

echo "deb http://packages.elastic.co/elasticsearch/1.5/debian stable main" | sudo tee -a /etc/apt/sources.list

sudo apt-get update && sudo apt-get install elasticsearch

sudo update-rc.d elasticsearch defaults 95 10

sudo service elasticsearch start


sudo apt-get install tcl8.5 -y
sudo apt-get install software-properties-common python-software-properties -y
sudo add-apt-repository ppa:chris-lea/redis-server -y
sudo apt-get update
sudo apt-get install redis-server -y




echo "create user chembl with password 'chembl'; create database cbh_chembl ; grant all privileges on database cbh_chembl to chembl;" > /tmp/creator

sudo su postgres -c "psql < /tmp/creator"
sudo su postgres -c "echo 'create extension if not exists hstore;create extension if not exists rdkit;'"
echo 'export DJANGO_SETTINGS_MODULE="deployment.settings/.development"' >> ~/.bashrc
echo 'source ~/anaconda/bin/activate chembiohub_wsa' >> ~/.bashrc
souce ~/.bashrc

  mkdir ~/.local/lib/python2.7/site-packages -p
  cd ~/.local/lib/python2.7/site-packages/
  wget https://raw.githubusercontent.com/thesgc/chembiohub_ws/master/scripts/vagrant_pythonpath.pth
pip install git+git://github.com/dcramer/django-devserver#egg=django-devserver

cd ~/chemreg/src/ng-chem/
sudo chmod o+w -R .
sudo npm install
bower install
