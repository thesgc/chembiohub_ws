cd /home/vagrant/chembiohub_ws/src/ng-chem
npm install
bower install
echo 'LOGIN_REDIRECT_URL="http://localhost:9000/#/projects/list" >> /home/vagrant/chembiohub_ws/deployment/settings/secret.py'