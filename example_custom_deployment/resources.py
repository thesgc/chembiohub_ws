import logging

# Get an instance of a logger
logger = logging.getLogger(__name__)

import subprocess
from pprint import pprint

from cbh_chembl_ws_extension.compounds import CBHCompoundBatchResource
from django.http import HttpRequest
from cbh_chembl_ws_extension.parser import get_uncurated_fields_from_file


def add_external_ids_to_file_object(python_file_obj):
    """Save the contents of the sdf file to the extenal system and fill in the external ID field
    Write the upadted file object to disc"""
    
    logger.debug("test log")
    logger.info("test 2 log")
    logger.error("test 3 log")

    logger.info(vars(python_file_obj))

    cmd_args = ['/usr/icm-3.8-4/icm64', '/home/chembiohub/scarab/ChemRegScarabPlugin.icm', '-a', 'sdf='+python_file_obj.path, 'runmode=GEN']

    subprocess.call(cmd_args)

    python_file_obj.path = python_file_obj.path + '_out'

    pass

def validate_file_object_externally(python_file_obj):
    """E.G. take the sdf file and set a status field to say whether it is in the external system
    Write the upadted file object to disc"""
    logger.info('Hello World')
    logger.info(vars(python_file_obj))

    cmd_args = ['/usr/icm-3.8-4/icm64', '/home/chembiohub/scarab/src/helpers/ChemRegHelper/ChemRegScarabPlugin.icm', '-a', 'sdf='+python_file_obj.name, 'runmode=GEN']

    subprocess.call(cmd_args)

    #python_file_obj.name = python_file_obj.name + '_out'

class AlteredCompoundBatchResource(CBHCompoundBatchResource):


    def alter_batch_data_after_save(self, batch_list, python_file_obj, request):
        """Take the original sdf file and run an external process with it such that new data can be written across 
        after save of the data into ChemBioHub"""
        cmd_args = ['/usr/icm-3.8-4/icm64', '/home/chembiohub/scarab/src/helpers/ChemRegHelper/ChemRegScarabPlugin.icm', '-a', 'sdf='+python_file_obj.name, 'runmode=INSERT']

        subprocess.call(cmd_args)


        add_external_ids_to_file_object(python_file_obj)
        fielderrors = {}
        fielderrors["stringdate"] = set([])
        fielderrors["number"] = set([])
        fielderrors["integer"] = set([])
        uncurated_data = get_uncurated_fields_from_file(python_file_obj, fielderrors)[0]
        for index, batch in enumerate(batch_list):
            #This assumes project is set up with exactly right custom fields
            batch.custom_fields = uncurated_data[index]
            batch.save()


    def preprocess_sdf_file(self, python_file_obj, request):
        """Hook to preprocess an sdf file - assumes that the file will be written to disc"""
        validate_file_object_externally(python_file_obj)
