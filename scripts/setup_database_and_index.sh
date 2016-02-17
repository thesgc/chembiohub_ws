
 #Start postgres in foreground


psql  -h $CONDA_ENV_PATH/var/postgressocket -c "create extension if not exists hstore;create extension if not exists  rdkit;" template1



createdb -h $CONDA_ENV_PATH/var/postgressocket/ ${ENV_NAME}_db -T template1

python manage.py migrate
python manage.py loaddata datatypes.json
python manage.py loaddata projecttypes.json
python manage.py reindex_compounds
python manage.py reindex_datapoint_classifications
python manage.py createsuperuser
python manage.py collectstatic
