
sudo useradd -G sudo -s /bin/bash -m chembl
sudo echo "chembl:chemblvm"| sudo chpasswd

sudo sysctl -w kernel.shmmax=2147483648
#sudo /sbin/iptables -A INPUT -i eth0 -p tcp --destination-port 5432 -j ACCEPT

#sudo dd if=/dev/zero of=/swapfile bs=1024 count=524288
#sudo chmod 600 /swapfile
#sudo mkswap /swapfile
#sudo swapon /swapfile

#wget https://raw.githubusercontent.com/chembl/mychembl/master/install_core_libs.sh && bash install_core_libs.sh
sudo apt-get install -y bash
sudo apt-get install -y git
sudo apt-get install -y unzip
#sudo apt-get install -y ipython
#sudo apt-get install -y ipython-notebook
#sudo apt-get install -y ipython-qtconsole
sudo apt-get install -y libboost-all-dev
sudo apt-get install -y postgresql
sudo apt-get install -y postgresql-server-dev-all
sudo apt-get install -y postgresql-doc
sudo apt-get install -y postgresql-contrib
sudo apt-get install -y flex
sudo apt-get install -y bison
sudo apt-get install -y g++
sudo apt-get install -y cmake
sudo apt-get install -y make
sudo apt-get install -y libffi-dev
sudo apt-get install -y libxml2 libxml2-dev
sudo apt-get install -y libxslt1.1 libxslt1-dev
sudo apt-get install -y python-numpy
sudo apt-get install -y python-scipy
sudo apt-get install -y python-matplotlib
sudo apt-get install -y python-pip
sudo apt-get install -y python-psycopg2
sudo apt-get install -y python-imaging-tk
sudo apt-get install -y python-networkx
sudo apt-get install -y libnss-mdns
sudo apt-get install -y avahi-utils
sudo apt-get install -y python-gobject
sudo apt-get install -y python-dev
sudo apt-get install -y python-biopython
sudo apt-get install -y rcconf



wget https://raw.githubusercontent.com/chembl/mychembl/master/install_py_libs.sh && bash install_py_libs.sh
wget https://raw.githubusercontent.com/chembl/mychembl/master/ensure_ipv6.sh && bash ensure_ipv6.sh

sudo -u postgres createuser -dsr chembl

cd /tmp
#sudo curl -o /etc/network/interfaces https://raw.githubusercontent.com/chembl/mychembl/master/configuration/mychembl_interfaces
#sudo curl -o /etc/init/failsafe.conf https://raw.githubusercontent.com/chembl/mychembl/master/configuration/failsafe.conf
sudo -u chembl curl -o /home/chembl/.bashrc https://raw.githubusercontent.com/chembl/mychembl/master/configuration/mychembl_bashrc
#sudo curl -o /etc/init/mychembl-upnp.conf https://raw.githubusercontent.com/chembl/mychembl/master/zeroconf/mychembl-upnp.conf
sudo curl -o /etc/avahi/services/mychembl.service https://raw.githubusercontent.com/chembl/mychembl/master/zeroconf/mychembl.service
sudo curl -o /usr/bin/mychembl-upnp.py https://raw.githubusercontent.com/chembl/mychembl/master/zeroconf/mychembl-upnp.py
sudo chmod +x /usr/bin/mychembl-upnp.py
sudo mkdir /usr/share/themes/mychembl
sudo curl -o /usr/share/themes/mychembl/mychembl.png https://raw.githubusercontent.com/chembl/mychembl/master/branding/mychembl.png
sudo curl -o /lib/plymouth/themes/ubuntu-text/ubuntu-text.plymouth https://github.com/chembl/mychembl/blob/master/branding/ubuntu-text.plymouth

cd /tmp
#yes
wget https://raw.githubusercontent.com/chembl/mychembl/master/osra.sh && bash osra.sh
wget https://raw.githubusercontent.com/chembl/mychembl/master/rdkit_install.sh && su -c "bash rdkit_install.sh" chembl

#wget https://raw.githubusercontent.com/chembl/mychembl/master/ipynb_setup.sh && su -c "bash ipynb_setup.sh" chembl
#Tidy up
#wget https://raw.githubusercontent.com/chembl/mychembl/master/create_db.sh && su -c "bash create_db.sh" chembl

#Self install
#wget https://raw.githubusercontent.com/chembl/mychembl/master/webservices/ws_setup.sh && su -c "bash ws_setup.sh" chembl
wget https://raw.githubusercontent.com/chembl/mychembl/master/beaker/install_beaker.sh && su -c "bash install_beaker.sh" chembl
#No
#wget https://raw.githubusercontent.com/chembl/mychembl/master/ipython_notebooks/ipynb_deamonise.sh && su -c "bash ipynb_deamonise.sh" chembl
#wget https://raw.githubusercontent.com/chembl/mychembl/master/apache.sh && bash apache.sh
#wget https://raw.githubusercontent.com/chembl/mychembl/master/launchpad.sh && bash launchpad.sh

sudo apt-get remove -y cloud-init
sudo rm /etc/init/cloud-* -rf

sudo swapoff -v /swapfile
sudo rm /swapfile
sudo rm /tmp/* -rf

#sudo bash -c "echo \"GRUB_BACKGROUND=\\\"/usr/share/themes/mychembl/mychembl.png\\\"\" >> /etc/default/grub"
sudo update-grub
