wget -N https://raw.githubusercontent.com/thesgc/chembiohub_ws/master/scripts/clone_project_to_given_folder.sh
sed -i -e '/Defaults\s\+env_reset/a Defaults\texempt_group=sudo' /etc/sudoers
sed -i -e 's/%sudo  ALL=(ALL:ALL) ALL/%sudo  ALL=NOPASSWD:ALL/g' /etc/sudoers

echo "UseDNS no" >> /etc/ssh/sshd_config
sudo -H -u vagrant bash -c "sh clone_project_to_given_folder.sh /srv/chembiohub vagrant"

