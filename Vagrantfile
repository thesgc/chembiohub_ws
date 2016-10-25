# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|

  config.vm.box = "https://chembiohub.ox.ac.uk/package.box"

  config.vm.network "forwarded_port", guest: 8000, host: 8000, auto_correct: true

  config.vm.network "forwarded_port", guest: 9000, host: 9000, auto_correct: true

  config.vm.network "forwarded_port", guest: 35729, host: 35729, auto_correct: true

  config.vm.provider "virtualbox" do |v|
    v.memory = 4096
    v.cpus = 2
    #v.customize(['storagectl', :id, '--name', 'SATAController', '--hostiocache', 'off'])
  end

  config.vm.synced_folder   ".", "/home/vagrant/syw" 
  #config.vm.provision :shell, :path => "" , :privileged => false
  config.ssh.password = "vagrant"
end
