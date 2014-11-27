chembiohub_ws
=============

A django project holder for the chembiohub web services project, submodules will be used for the app projects which are in concurrent development.

The concept of the project is to extend chembl_core_model and chembl_webservices to provide additional functionality. This will then be exposed in a separate angularjs app.

The web service and model component repositories were added as subrepositories to allow easy development but at the same time keep the code foir web services separate to that for models

This was done by first running:

  git submodule add git@github.com:thesgc/cbh_chembl_model_extension
  git submodule add git@github.com:thesgc/cbh_chembl_ws_extension

As the modules are full coookiecutter installable apps then the django app folder is one level down from the initial directory, hence in settings.py the installed app is for example:

  cbh_chembl_model_extension.cbh_chembl_model_extension

With a small change to settings.py these will eventually be installable using pip.

Development is being done in the chembl 19 virtual machine using the existing virtual environement to begin with.
