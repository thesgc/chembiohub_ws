chembiohub_ws
=============

The main holder project for ChemBio Hub Platform including back end web services, database and front end code all as separate submodules

All relevant code is in the /src/ directory as sub repositories

Current install instructions do not work on Windows but work on Mac and Linux.

In order to get started using vagrant then run to get all the subrepos:

    git clone  --recursive  git@github.com:thesgc/chembiohub_ws.git
  
Next we create a vagrant machine to run our python code:

    vagrant up
  
This will download the vagrant box from our internet location and you can then log in by:

    vagrant ssh
  
You should now see the prompt to show that the  virtualenv is enabled and you are logged on to the vagrant box. 


To get the front end dependencies run:

    cd ~/chembiohub_ws/src/ng-chem
    sudo npm install && bower install

In order to migrate the database to the latest version run:

    cd ~/chembiohub_ws/
    python manage.py migrate && python manage.py collectstatic

To run the back end server run

    cd ~/chembiohub_ws/
    python manage.py runserver 0.0.0.0:8000

You will now have the server running inside the vagrant box. 

Next, open another terminal window and run

    vagrant ssh
    cd  ~/chembiohub_ws/src/ng-chem
    grunt serve

You can then access the login page via

    http://localhost:9000/dev/login

The username and password are both

    vagrant


For a production installation please get in touch and we will help you get apache or nginx configured, guides are coming soon.



Alternatively to run the grunt on your local box open a second terminal window and change directory to the ng-chem repository

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




All of these are in development mode only, for a production install you will need to set up an appropriate settings file and follow the instructions in [install ubuntu](install_ubuntu.md) or [install centos](install_centos.md).


