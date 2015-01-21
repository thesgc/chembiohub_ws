
wget http://www.iupac.org/fileadmin/user_upload/publications/e-resources/inchi/1.03/INCHI-1-BIN.zip
unzip INCHI-1-BIN.zip
gunzip INCHI-1-BIN/linux/64bit/inchi-1.gz

cd chembiohub_ws
wget https://dl.dropboxusercontent.com/u/10967207/indigo-python-1.1.11-linux.zip

unzip indigo-python-1.1.11-linux.zip

rm indigo-python-1.1.11-linux.zip



echo "" >> ~/.bashrc 
echo 'export DJANGO_SETTINGS_MODULE="deployment.settings.staging"' >> ~/.bashrc

source ~/miniconda/bin/activate newbeak
pip install django-cors-headers
cd ~/chembiohub_ws/src/chembl_core_db
python setup.py install
python setup.py develop

cd ../chembl_core_model/
python setup.py install
python setup.py develop

cd ../chembl_webservices/
python setup.py install
python setup.py develop

cd ../chembl_business_model/
python setup.py install
python setup.py develop

cd ../standardiser/
python setup.py install
python setup.py develop

cd ../chembl_beaker/
python setup.py install
python setup.py develop

cd ../cbh_chembl_model_extension/
python setup.py install
python setup.py develop

cd ../cbh_chembl_ws_extension/
python setup.py install
python setup.py develop

cd ../chembl_extras/
python setup.py install
python setup.py develop

sudo apt-get install nodejs -y

sudo apt-get install npm -y

sudo npm install -g bower

sudo apt-get install nodejs-legacy -y

cd /home/vagrant/chembiohub_ws/src/ng-chem

bower install