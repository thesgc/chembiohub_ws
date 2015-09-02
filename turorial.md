#How to configure ChemBio Hub ChemReg

##Local installation:

First install vagrant and virtualbox locally

Install the codebase using vagrant like this:

    git clone https://github.com/thesgc/chembiohub_ws
    git submodule init
    git submodule update
    vagrant up
    
Once the vagrant machine has deployed and started, you will need to connect in 2 windows in order to run ChemBio Hub ChemReg

Open 2 terminal windows to the source folder and run:

    vagrant ssh

In each terminal

In one termiinal we will start the python server, in the other, the javascript server:

    cd ~/chembiohub_ws python manage.py runserver

In a second terminal run

    cd ~/chembiohub_ws/src/ng-chem
    grunt serve

You should now be able to log on to the server and set up ChemBio Hub ChemReg

First log on to the system via:

http://localhost:8000/dev/login

You will see a list of 3 projects

For this tutorial we are going to import some zinc data and add a plugin to add a custom descriptor to our data.

