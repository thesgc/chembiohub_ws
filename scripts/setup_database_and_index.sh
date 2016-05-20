
 #Start postgres in foreground
export ENV_NAME="$1"
OLD_PATH="$PATH"
mkdir -p /home/$2/anaconda2/envs/$ENV_NAME/var/postgressocket/

sudo supervisorctl reread
sudo supervisorctl update
sudo supervisorctl reload

printf "sleeping"

sleep 10
printf "finished sleeping"
#source activate $ENV_NAME

export CONDA_ENV_PATH=/home/$2/anaconda2/envs/$ENV_NAME
export PATH=$CONDA_ENV_PATH/bin:$OLD_PATH


$CONDA_ENV_PATH/bin/psql  -h $CONDA_ENV_PATH/var/postgressocket -c "create extension if not exists hstore;create extension if not exists  rdkit;" template1
$CONDA_ENV_PATH/bin/createdb -h $CONDA_ENV_PATH/var/postgressocket ${ENV_NAME}_db -T template1

$CONDA_ENV_PATH/bin/python deployment/generate_secret_settings.py > deployment/settings/secret.py
$CONDA_ENV_PATH/bin/python manage.py migrate
$CONDA_ENV_PATH/bin/python manage.py loaddata projecttypes.json
$CONDA_ENV_PATH/bin/python manage.py reindex_compounds_new
cd src/ng-chem
bower install
cd ../..
if [ "$ENV_NAME" != "dev" ]
then
    $CONDA_ENV_PATH/bin/python manage.py collectstatic --noinput
fi
echo "from django.contrib.auth.models import User; User.objects.create_superuser('admin', 'admin@example.com', 'pass')" | $CONDA_ENV_PATH/bin/python manage.py shell
