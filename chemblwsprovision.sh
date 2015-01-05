cd chembiohub_ws
wget https://dl.dropboxusercontent.com/u/10967207/indigo-python-1.1.11-linux.zip

unzip indigo-python-1.1.11-linux.zip

rm indigo-python-1.1.11-linux.zip

echo "export PYTHONPATH=:/var/www/chembiohub_ws:/home/chembiohub/indigo-python-1.1.11-linux:/home/chembiohub/Tools/openbabel-install/lib" >> ~/.bashrc 
echo 'export DJANGO_SETTINGS_MODULE="deployment.settings.staging"' >> ~/.bashrc

source ~/miniconda/bin/activate newbeak
pip install django-cors-headers
cd ~/chembiohub_ws/src/chembl_core_db
python setup.py develop

cd ../chembl_core_model/

python setup.py develop

cd ../chembl_webservices/

python setup.py develop

cd ../chembl_business_model/

python setup.py develop

cd ../standardiser/

python setup.py develop

cd ../chembl_beaker/

python setup.py develop

cd ../cbh_chembl_model_extension/

python setup.py develop

cd ../cbh_chembl_ws_extension/

python setup.py develop

cd ../chembl_extras/

python setup.py develop