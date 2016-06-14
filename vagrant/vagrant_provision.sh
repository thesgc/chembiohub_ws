#ln -s /srv/chembiohub /home/vagrant/syw/code
cd /srv/chembiohub
bash /srv/chembiohub/scripts/setup_database_and_index.sh dev vagrant
echo "LOGIN_REDIRECT_URL = '/#/projects/list'" >> /srv/chembiohub/deployment/settings/secret.py
