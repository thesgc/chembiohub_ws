cd /srv/chembiohub/
echo "update complete"
sudo -H -u vagrant bash -c "bash scripts/setup_database_and_index.sh chembiohub vagrant"