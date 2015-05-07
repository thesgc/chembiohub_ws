chembiohub_ws
=============

A django project holder for the chembiohub web services project, submodules will be used for the app projects which are in concurrent development.

The concept of the project is to extend chembl_core_model and chembl_webservices to provide additional functionality. This will then be exposed in a separate angularjs app.

All relevant code is in the /src/ directory as sub repositories

These are installed in your local anaconda install as shown in  [install anaconda](install_anaconda.rst) but this has been done for you on the vagrant install.

In order to get started using vagrant then run to get all the subrepos:

    git clone  --recursive  git@github.com:thesgc/chembiohub_ws.git
  
Next we create a vagrant machine to run our python code:

    vagrant up
  
This will download the vagrant box from our internet location and you can then log in by:

    vagrant ssh
  
You should now see the prompt to show that the  virtualenv is enabled. Addtionally, openbabel python libraries are on the classpath.

    cd ~/chembiohub_ws/
    python manage.py runserver 0.0.0.0:8000

You will now have the server running inside the vagrant box.

In order to take advantage of live reload on the front end then we use grunt serve for development.

On your local box open a second terminal window and change directory to the ng-chem repository
    cd src/ng-chem

Install the bower dependencies using the following for an ubuntu machine
   sudo apt-get install -y nodejs
  sudo apt-get install -y npm
  sudo apt-get install -y nodejs-legacy
  sudo apt-get install -y ruby gem ruby-dev
  sudo gem install compass

  sudo npm install -g bower grunt-cli coffee-script

You can then run
   npm install
   bower install

You can then run the server using

   grunt serve
   
This will allow the server to run locally with live reload on port 9000

In order to create a superuser run:

   python manage.py createsuperuser
   
in the vagrant propmt

Log in to the site by going to the login URL at:

    http://localhost:9000/dev/login

Add a project for your new user by going to the admin URL at
    localhost:8000/dev/admin
Next add a custom field config for the project

    
Next go back to localhost:9000/devapi/login and you will be redirected to the projects list with the new project in it

As creator of the project you will have permissions for the project.

Other project permissions can be edited on a per-user or per group basis




All of these are in development mode only, for a production install you will need to set up an appropriate settings file and follow the instructions in [install anaconda](install_anaconda.rst) 


