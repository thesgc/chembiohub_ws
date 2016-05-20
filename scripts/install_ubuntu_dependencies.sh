sudo apt-get install -y openjdk-7-jre-headless
wget -qO - https://packages.elastic.co/GPG-KEY-elasticsearch | sudo apt-key add -
echo "deb http://packages.elastic.co/elasticsearch/2.x/debian stable main" | sudo tee -a /etc/apt/sources.list.d/elasticsearch-2.x.list
sudo add-apt-repository ppa:chris-lea/node.js -y
sudo apt-get update && sudo apt-get install elasticsearch nodejs
sudo update-rc.d elasticsearch defaults 95 10
sudo service elasticsearch start
sudo add-apt-repository ppa:chris-lea/redis-server -y
sudo apt-get update
sudo apt-get install -y redis-server  ruby gem ruby-dev unzip fabric git apache2 libapache2-mod-proxy-html  libxml2-dev supervisor  tcl8.5 software-properties-common python-software-properties 
sudo update-rc.d redis-server defaults
sudo service redis-server start

if [ "$USER" != "travis" ]
    then
    sudo a2enmod proxy proxy_http headers expires rewrite
    sudo gem install compass
    sudo npm install -g bower grunt-cli coffee-script
fi
