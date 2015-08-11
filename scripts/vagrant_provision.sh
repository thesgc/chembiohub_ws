echo "create user chembl with password 'chembl'; create database cbh_chembl ; grant all privileges on database cbh_chembl to chembl;" > /tmp/creator

sudo su postgres -c "psql < /tmp/creator"
sudo su postgres -c "echo 'create extension if not exists hstore;create extension if not exists rdkit;'"
echo 'export DJANGO_SETTINGS_MODULE="deployment.settings.development"' > ~/.bashrc

source ~/.bashrc

  mkdir ~/.local/lib/python2.7/site-packages -p
  cd ~/.local/lib/python2.7/site-packages/
  wget https://raw.githubusercontent.com/thesgc/chembiohub_ws/master/scripts/vagrant_pythonpath.pth
pip install django_devserver
