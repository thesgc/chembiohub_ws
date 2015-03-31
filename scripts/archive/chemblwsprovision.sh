
wget http://www.iupac.org/fileadmin/user_upload/publications/e-resources/inchi/1.03/INCHI-1-BIN.zip
unzip INCHI-1-BIN.zip
gunzip INCHI-1-BIN/linux/64bit/inchi-1.gz

cd chembiohub_ws
wget https://dl.dropboxusercontent.com/u/10967207/indigo-python-1.1.11-linux.zip

unzip indigo-python-1.1.11-linux.zip

rm indigo-python-1.1.11-linux.zip


echo 'source ~/miniconda/bin/activate chembiohub_ws' >> /home/vagrant/.bashrc
echo 'export PYTHONPATH=:/home/vagrant/chembiohub_ws:/home/vagrant/chembiohub_ws/indigo-python-1.1.11-linux:/home/vagrant/Tools/openbabel-install/lib:/home/vagrant/chembiohub_ws:/home/vagrant/chembiohub_ws/src/chembl_core_db/:/home/vagrant/chembiohub_ws/src/chembl_core_model/:/home/vagrant/chembiohub_ws/src/chembl_business_model/:/home/vagrant/chembiohub_ws/src/standardiser/:/home/vagrant/chembiohub_ws/src/chembl_beaker/:/home/vagrant/chembiohub_ws/src/cbh_chembl_model_extension/:/home/vagrant/chembiohub_ws/src/cbh_chembl_ws_extension/:/home/vagrant/chembiohub_ws/src/django-flowjs/' >> /home/vagrant/.bashrc
echo 'export DJANGO_SETTINGS_MODULE="deployment.settings.development"' >> /home/vagrant/.bashrc

