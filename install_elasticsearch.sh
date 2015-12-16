wget https://download.elasticsearch.org/elasticsearch/release/org/elasticsearch/distribution/deb/elasticsearch/2.1.0/elasticsearch-2.1.0.deb
sudo dpkg -i elasticsearch-2.1.0.deb
sudo service elasticsearch start
sudo update-rc.d elasticsearch defaults
sudo /usr/share/elasticsearch/bin/elasticsearch install license
sudo /usr/share/elasticsearch/bin/elasticsearch install shield
