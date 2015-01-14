from __future__ import with_statement
from fabric.api import local
from fabric.api import *
from fabric.contrib.console import confirm
import os
env.hosts = ['staging.chembiohub.ox.ac.uk']

def _deploy(code_dir, process_name):
    # with settings(warn_only=True):
    #     if run("test -d %s" % code_dir).failed:
    #         run("git clone user@vcshost:/path/to/repo/.git %s" % code_dir)
    with cd(code_dir):          
        sudo("git pull origin master", user="chembiohub") 
        sudo("git submodule init", user="chembiohub") 
        sudo("git submodule update", user="chembiohub") 
        sudo("git submodule foreach git pull origin master", user="chembiohub")        
        sudo("supervisorctl reload")
        sudo("service apache2 reload reload")
        with cd("src/ng-chem"):
            run("bower install") 

        for dirname in dirs[1:]:
            with cd("src/%s" % dirname):
                sudo("source /home/chembiohub/.bashrc &&  source /home/chembiohub/anaconda/bin/activate chembiohub_ws && python setup.py install ", user="chembiohub")
        sudo("source /home/chembiohub/.bashrc &&  source /home/chembiohub/anaconda/bin/activate chembiohub_ws && pip install -r requirements_additional.txt && python manage.py migrate cbh_chembl_model_extension &&   python manage.py syncdb  && python manage.py collectstatic", user="chembiohub") 



def prod():
    # env.hosts = ['chembiohub.ox.ac.uk']
    _deploy("/var/www/chembiocrunch/", "chem_bio_crunch")


def stage():
    
    _deploy("/var/www/chembiohub_ws", "chem_bio_hub_ws")


dirs = ["ng-chem",
 "chembl_core_db",
        "chembl_core_model",
 "chembl_business_model",
        
        "chembl_extras",
        "chembl_webservices",
        "cbh_chembl_model_extension",
                "django-flowjs",
        "standardiser",
        "chembl_beaker",

        "cbh_chembl_ws_extension",
       
        
        ]


def prep():
    count = 4
    _prep(count)

def prepall():
    _prep(len(dirs))


def _prep(dircount):
    dirnow = os.getcwd()
    for directory in dirs[0:dircount-1]:

        with lcd(dirnow + "/src/" + directory):
            print(directory)
            try:
                local("git add  .")
            except:
                pass
            try:
                local("git commit -a")
            except:
                pass
            try:
                local("git pull origin master")
            except:
                pass
            try: 
                local("git push origin master")
            except:
                pass
    try:
        local("git add .")
    except:
        pass
    try:
        local("git commit -a")
    except:
        pass
    try:
        local("git pull origin master")
    except:
        pass
    try: 
        local("git push origin master")
    except:
        pass

