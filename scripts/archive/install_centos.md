
First install redis, supervisor and elasticsearch and apache2

    yum install wget gcc python-meld3
    wget http://download.redis.io/releases/redis-2.8.3.tar.gz
    tar xzvf redis-2.8.3.tar.gz
    cd redis-2.8.3
    make
    make install
    rpm -Uvh http://download.fedoraproject.org/pub/epel/6/x86_64/epel-release-6-8.noarch.rpm
    rpm -Uvh http://rpms.famillecollet.com/enterprise/remi-release-6.rpm
    yum --enablerepo=remi,remi-test install redis
    echo "sysctl vm.overcommit_memory=1" >> /etc/sysctl.conf
    sysctl -w fs.file-max=100000
    chkconfig --add redis
    chkconfig --level 345 redis on
    service redis start/stop/restart
    service redis start

    sudo yum install python-setuptools
    sudo easy_install pip
    sudo pip install supervisor
    mkdir -p /etc/supervisor/conf.d
    echo_supervisord_conf > supervisord.conf
    sudo mv supervisord.conf /etc/supervisor/supervisord.conf
    printf "[include]\nfiles  =  /etc/supervisor/conf.d/*.conf" >> /etc/supervisor/supervisord.conf
    printf '#!/bin/sh\n#\n# /etc/rc.d/init.d/supervisord\n#\n# Supervisor is a client/server system that\n# allows its users to monitor and control a\n# number of processes on UNIX-like operating\n# systems.\n#\n# chkconfig: - 64 36\n# description: Supervisor Server\n# processname: supervisord\n\n# Source init functions\n. /etc/rc.d/init.d/functions\n\nprog="supervisord"\n\nprefix="/usr/"\nexec_prefix="${prefix}"\nprog_bin="${exec_prefix}/bin/supervisord -c /etc/supervisor/supervisord.conf"\nPIDFILE="/var/run/$prog.pid"\n\nstart()\n{\n       echo -n $"Starting $prog: "\n       daemon $prog_bin --pidfile $PIDFILE\n       [ -f $PIDFILE ] && success $"$prog startup" || failure $"$prog startup"\n       echo\n}\n\nstop()\n{\n       echo -n $"Shutting down $prog: "\n       [ -f $PIDFILE ] && killproc $prog || success $"$prog shutdown"\n       echo\n}\n\ncase "$1" in\n\n start)\n   start\n ;;\n\n stop)\n   stop\n ;;\n\n status)\n       status $prog\n ;;\n\n restart)\n   stop\n   start\n ;;\n\n *)\n   echo "Usage: $0 {start|stop|restart|status}"\n ;;\n\nesac' > /etc/rc.d/init.d/supervisord
    sudo chmod +x /etc/rc.d/init.d/supervisord
    sudo chkconfig --add supervisord
    sudo chkconfig supervisord on
    sudo service supervisord start



    cd ~
    wget --no-cookies --no-check-certificate --header "Cookie: gpw_e24=http%3A%2F%2Fwww.oracle.com%2F; oraclelicense=accept-securebackup-cookie" "http://download.oracle.com/otn-pub/java/jdk/8u65-b17/jdk-8u65-linux-x64.rpm"
    sudo yum localinstall jdk-8u65-linux-x64.rpm
    rm ~/jdk-8u65-linux-x64.rpm
    sudo rpm --import http://packages.elastic.co/GPG-KEY-elasticsearch
    printf "[elasticsearch-2.0]\nname=Elasticsearch repository for 2.x packages\nbaseurl=http://packages.elastic.co/elasticsearch/2.x/centos\ngpgcheck=1\ngpgkey=http://packages.elastic.co/GPG-KEY-elasticsearch\nenabled=1\n\n\n" > /etc/yum.repos.d/elasticsearch.repo
    sudo yum -y install elasticsearch
    sudo systemctl start elasticsearch
    chkconfig --add elasticsearch
    chkconfig --level 345 elasticsearch on
    service elasticsearch start


    yum install httpd
    sudo service httpd start
    sudo chkconfig httpd on
    iptables -I INPUT -p tcp -m tcp --dport 80 -j ACCEPT
    /sbin/service iptables save

    printf "# This file controls the state of SELinux on the system.
# SELINUX= can take one of these three values:
#     enforcing - SELinux security policy is enforced.
#     permissive - SELinux prints warnings instead of enforcing.
#     disabled - No SELinux policy is loaded.
SELINUX=disabled
# SELINUXTYPE= can take one of these two values:
#     targeted - Targeted processes are protected,
#     mls - Multi Level Security protection.
SELINUXTYPE=targeted" > /etc/selinux


Next create a chembiohub user and switch to it

    sudo adduser chembiohub
    sudo su chembiohub
    cd ~

    git clone https://github.com/thesgc/chembiohub_ws
    cd chembiohub_ws
    git checkout conda
    bash scripts/install_linux64.sh 0 Centos # Zero is the directory this will be installed into
