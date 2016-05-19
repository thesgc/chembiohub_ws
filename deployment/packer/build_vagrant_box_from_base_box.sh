cd /srv/chembiohub
git pull
git submodule foreach git pull origin master
bash -c "bash scripts/setup_database_and_index.sh dev vagrant"
bash -c "bash scripts/run_bower.sh dev Ubuntu"



