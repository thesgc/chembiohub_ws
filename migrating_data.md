#In order to mogitate data from one instance of chembiohub to another the database and potentially index can be backed up and migrated.

To migrate the database we use pg_dump and pg_restore
As postgres user
    pg_dump -Fc hub_reg_db > /tmp/hub_reg.bak

As chembiohub user
    pg_restore -Fc -d hubtest_db /tmp/hub_reg.bak -h $CONDA_ENV_PATH/var/postgressocket

To migrate the index we use the elasticdump tool

elasticdump \
  --input=http://localhost:9200/hub__cbh_datastore_index \
  --output=http://localhost:9200/hubtest__cbh_datastore_index \
  --type=mapping
elasticdump \
   --input=http://localhost:9200/hub__cbh_datastore_index \
  --output=http://localhost:9200/hubtest__cbh_datastore_index \
  --type=data


elasticdump \
  --input=http://localhost:9200/hub__chemreg_chemical_index \
  --output=http://localhost:9200/hub__chemreg_chemical_index \
  --type=mapping
elasticdump \
  --input=http://localhost:9200/hub__chemreg_chemical_index \
  --output=http://localhost:9200/hub__chemreg_chemical_index \
  --type=data
