from setuptools import setup, find_packages
setup(
    name = "ShowYourWorking",
    version = "0.1",
    packages = ["deployment","cbh_core_model", "cbh_chembl_model_extension", 
   'chembl_core_db',
   'chembl_core_model',
  'chembl_business_model', 
         'cbh_core_api',
       'cbh_core_model',
  'cbh_chembl_model_extension',    
  'cbh_chem_api',
        'cbh_tests',
'cbh_utils',
          ],
  #  scripts = ['say_hello.py'],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
 #   install_requires = ['docutils>=0.3'],

    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
        # And include any *.msg files found in the 'hello' package, too:
        'hello': ['*.msg'],
    },

    # metadata for upload to PyPI
    author = "Me",
    author_email = "me@example.com",
    description = "This is an Example Package",
    license = "PSF",
    keywords = "hello world example examples",
    url = "http://example.com/HelloWorld/",   # project home page, if any

    # could also include long_description, download_url, classifiers, etc.
)
