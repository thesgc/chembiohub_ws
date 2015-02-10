prefix = "test"
port = "8081"
base_folder = "/var/www"


apache_template = '''   ProxyPass /{prefix}/reg !
    Alias /{prefix}/reg {base_folder}/{prefix}/chembiohub_ws/deployment/static/dist

    <Directory {base_folder}/{prefix}/chembiohub_ws/deployment/static/dist>
        Options Indexes FollowSymLinks Includes
        AllowOverride All
        Order allow,deny
        Allow from all
    </Directory>
  ProxyPass /{prefix}/ http://localhost:{port}/{prefix}/
  ProxyPassReverse /{prefix}/ http://localhost:{port}/{prefix}/
'''
