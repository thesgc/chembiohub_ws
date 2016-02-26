set -e
 #Start postgres in foreground
export ENV_NAME="$1"
OLD_PATH="$PATH"


#source activate $ENV_NAME

export CONDA_ENV_PATH=/home/$2/anaconda2/envs/$ENV_NAME
export PATH=$CONDA_ENV_PATH/bin:$OLD_PATH


$CONDA_ENV_PATH/bin/psql  -h $CONDA_ENV_PATH/var/postgressocket -c "create extension if not exists hstore;create extension if not exists  rdkit;" template1
$CONDA_ENV_PATH/bin/createdb -h $CONDA_ENV_PATH/var/postgressocket/ ${ENV_NAME}_db -T template1

$CONDA_ENV_PATH/bin/python generate_secret_settings.py > deployment/settings/secret.py
$CONDA_ENV_PATH/bin/python manage.py migrate
$CONDA_ENV_PATH/bin/python manage.py loaddata datatypes.json
$CONDA_ENV_PATH/bin/python manage.py loaddata projecttypes.json
$CONDA_ENV_PATH/bin/python manage.py reindex_compounds
$CONDA_ENV_PATH/bin/python manage.py reindex_datapoint_classifications
cd src/ng-chem
bower install
cd ../..
$CONDA_ENV_PATH/bin/python manage.py collectstatic --noinput
echo "from django.contrib.auth.models import User; User.objects.create_superuser('admin', 'admin@example.com', 'pass')" | $CONDA_ENV_PATH/bin/python manage.py shell
