import logging

# Get an instance of a logger
logger = logging.getLogger(__name__)

import subprocess
from pprint import pprint

from cbh_chem_api.compounds import CBHCompoundUploadResource
from cbh_chem_api.resources import BaseCBHCompoundBatchResource
from django.http import HttpRequest
from cbh_utils.parser import get_uncurated_fields_from_file
from copy import copy
import json

from tastypie.exceptions import BadRequest

from datetime import date

import json
import os
import pprint

import traceback
import re

def add_external_ids_to_file_object(python_file_obj):
	"""Save the contents of the sdf file to the extenal system and fill in the external ID field
	Write the upadted file object to disc"""

def validate_file_object_externally(python_file_obj):
	"""E.G. take the sdf file and set a status field to say whether it is in the external system
	Write the upadted file object to disc"""
	logger.info('File ' + python_file_obj.name)
   
	cmd_args = ['/usr/icm-3.8-4/icm64', '/home/chembiohub/scarab/src/helpers/ChemRegHelper/ChemRegScarabPlugin.icm', '-a', 'sdf='+python_file_obj.name, 'runmode=GEN', 'username=AUTO']

	logger.info('File ' + python_file_obj.name)

	try:
		process = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		output = process.communicate()
	except error:
	logger.info(output[0])

		logger.info('An error ocurred running the registration code' + error)
		raise BadRequest('An unexpected error has occurred')

	m = re.search('<Error>([^<]+)',output[0])
	if not m is None:
		error = m.group(1)
	logger.info('Raising exception ' + error)
		raise BadRequest('An error has occurred ' + error)
	else:
	m = re.search('<Success>', output[0])
	logger.info(output[0])

	if m is None:
		logger.info('Something has gone wrong registering compounds')
		raise BadRequest('An unexpected error has occurred')

class AlteredCompoundBatchResource(CBHCompoundUploadResource):
	def after_save_and_index_hook(self, request, multi_batch_id, project_id):
		logger.info('Writing out files to share')
		try:
            extra_queries = [{'query_type': 'pick_from_list', 'field_path': 'multiple_batch_id','pick_from_list': [str(multiple_batch_id)]}]
			logger.info('Running Excel fetch')
			# Get Excel 
			newrequest = copy(request)
			newrequest.GET = request.GET.copy()
			newrequest.GET["format"] = "xlsx"
			newrequest.GET["multiple_batch_id"] = multi_batch_id
			newrequest.GET["offset"] = 0
			newrequest.GET["limit"] = 10000
			newrequest.GET["pids"] = str(project_id)
			file_resp = BaseCBHCompoundBatchResource().get_list(newrequest2, extra_queries=extra_queries)

			# Fetch username
			username = request.user.username

			logger.info('Running JSON fetch')
			# Fetch data to get prefix
			newrequest2 = copy(request)
			newrequest2.GET = request.GET.copy()
			newrequest2.GET['FORMAT'] = 'json'
			newrequest2.GET["multiple_batch_id"] = multi_batch_id
			newrequest2.GET["offset"] = 0
			newrequest2.GET["limit"] = 10000
			newrequest2.GET["pids"] = str(project_id)
			file_resp = BaseCBHCompoundBatchResource().get_list(newrequest2, extra_queries=extra_queries)
			ser_json = resp._container[0]
			json_object = json.loads(ser_json)

			global_prefix = json_object["objects"][0]["customFields"]["RequestedGlobalId"]

			# Fetch prefix		

			logger.info('Building Excel path')
			# Write out Excel
			excel_path = '/home/share/ChemInformatics/Registrations/' + username + '_' + global_prefix + '_'  + date.today().strftime('%d%m%y') + '.xlsx'

			# Handle multiple registration requests from the same user on the same day for the same prefix
			i = 1

			logger.info('Working out Excel file name')
			while(True):
				if(not os.path.exists(excel_path)):
				break
				else:
				i += 1

				if(i == 2):
					excel_path = excel_path + '_'
			
				excel_path = excel_path + str(i)

			logger.info('Writing out ' + excel_path)
			with open(excel_path, 'wb') as f:
				f.write(file_resp._container[0])
		 
			logger.info('Fetching SDF')
			# Fetch as SDF
			newrequest = copy(request)
			newrequest.GET = request.GET.copy()
			newrequest.GET["format"] = "sdf"
			newrequest.GET["multiple_batch_id"] = multi_batch_id
			newrequest.GET["offset"] = 0
			newrequest.GET["limit"] = 10000
			newrequest.GET["pids"] = str(project_id)
			file_resp = BaseCBHCompoundBatchResource().get_list(newrequest, extra_queries=extra_queries)

			# Output SDF
			sdf_path = '/home/share/ChemInformatics/Registrations/' + username + '_' + global_prefix + '_'  + date.today().strftime('%d%m%y') + '.sdf'

			# Handle multiple registrations requests from the same user on the same day for the same prefix
			i = 1
			
			logger.info('Working out SDF file name')
			while(True):
				if(not os.path.exists(sdf_path)):
				break
				else:
				i += 1

				if(i == 2):
					sdf_path = sdf_path + '_'
			
				sdf_path = sdf_path + str(i)

			logger.info('Writing out ' + sdf_path)
			with open(sdf_path, 'wb') as f:
				f.write(file_resp._container[0])
		except Exception as e:
			logger.info(traceback.format_exc())
			logger.info('An error has occurred' + e.__doc__ + '/' + e.message)
		logger.info('Export done')

	def alter_batch_data_after_save(self, batch_list, python_file_obj, request, multi_batch):
		"""Take the original sdf file and run an external process with it such that new data can be written across 
		after save of the data into ChemBioHub"""

		cmd_args = ['/usr/icm-3.8-4/icm64', '/home/chembiohub/scarab/src/helpers/ChemRegHelper/ChemRegScarabPlugin.icm', '-a', 'sdf='+python_file_obj.name, 'runmode=INSERT', 'username=AUTO']
		try:
			process = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			output = process.communicate()
		except error:
			logger.info('An error ocurred running the registration code' + error)
			raise BadRequest('An unexpected error has occurred')

		m = re.search('<Error>([^<]+)',output[0])
		if not m is None:
		logger.info(output[0])
			error = m.group(1)
		logger.info('Raising exception ' + error)
			raise BadRequest('An error has occurred ' + error)
		else:
		m = re.search('<Success>', output[0])
			logger.info(output[0])

			if m is None:
			logger.info('Something has gone wrong registering compounds')
			raise BadRequest('An unexpected error has occurred')

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

	def preprocess_sdf_file(self, python_file_obj,request):
		"""Hook to preprocess an sdf file - assumes that the file will be written to disc"""
		validate_file_object_externally(python_file_obj)
