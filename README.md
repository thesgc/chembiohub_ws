chembiohub_ws
=============

A django project holder for the chembiohub web services project, submodules will be used for the app projects which are in concurrent development.

The concept of the project is to extend chembl_core_model and chembl_webservices to provide additional functionality. This will then be exposed in a separate angularjs app.

All relevant code is in the /src/ directory as sub repositories

These are installed in your local anaconda install as shown in  [install anaconda](install_anaconda.rst) but this has been done for you on the vagrant install.

In order to get started using vagrant then run to get all the subrepos:

    git clone  --recursive  git@github.com:thesgc/chembiohub_ws.git
  
Next run:

    vagrant up
  
This will download the vagrant box from our internet location and you can then log in by:

    vagrant ssh
  
You should now see the prompt to show that the newbeak virtualenv is enabled. Addtionally, openbabel python libraries are on the classpath.

It should now be possible to run

    run_beaker
  
In order to see the beaker web services

Addtionally, to see the chembiohub version of the chembl web services, you can run

    cd ~/chembiohub_ws/deployment/
    python manage.py runserver
  
All of these are in development mode only, for a production install you will need to set up an appropriate settings file and follow the instructions in [install anaconda](install_anaconda.rst) 


