
 #Start postgres in foreground
export ENV_NAME="$1"
OLD_PATH="$PATH"


#source activate $ENV_NAME

export CONDA_ENV_PATH=$(conda info | grep "envs dir" | cut -c 25-)/$ENV_NAME
export PATH=$CONDA_ENV_PATH/bin:$OLD_PATH


createdb -h $CONDA_ENV_PATH/var/postgressocket/ ${ENV_NAME}_db -T template1
python generate_secret_settings.py > deployment/settings/secret.py
python manage.py migrate
python manage.py loaddata datatypes.json
python manage.py loaddata projecttypes.json
python manage.py reindex_compounds
python manage.py reindex_datapoint_classifications
cd src/ng-chem
bower install
cd ../..
python manage.py collectstatic --noinput
echo "from django.contrib.auth.models import User; User.objects.create_superuser('admin', 'admin@example.com', 'pass')" | python manage.py shell
