# -*- coding: utf-8 -*-
"""This module provides extra serializers to convert data into SDF or XLSX from the new search api
unfinished refactoring"""
from StringIO import StringIO

from tastypie.serializers import Serializer

import re
import json
import pandas as pd


import base64
import copy
import os

from cbh_core_api.resources import get_field_name_from_key
from cbh_core_api.resources import get_key_from_field_name
from cbh_core_api import parser as coreparser

from tastypie.exceptions import ImmediateHttpResponse, BadRequest
from jsonpointer import resolve_pointer
import xlsxwriter

SDF_TEMPLATE = u">  <{name}>\n{value}\n\n"

EMPTY_FILE_DESCRIPTION = "No files"
# Adding hyperlinks is not really possible because the Excel format only really
ATTACHMENT_TEMPLATE = ' {base_url}{url} {printName}'
#EXCEL_FORMAT  = '=HYPERLINK("{base_url}/#/searchv2/{doc_id}", "{number} files Ctrl + click to view/download")'

EXCEL_FORMAT  = '{base_url}/#/searchv2/{doc_id}'

def prepare_file_output(input_data, schema, export_type, doc_id):

    if isinstance(input_data, dict):
        list_of_attachments = input_data.get("attachments", [])
        if len(list_of_attachments) > 0:
            if export_type == "xlsx":
                return EXCEL_FORMAT.format(**{"base_url": schema.get("base_url", ""),
                                              "number" : len(list_of_attachments),
                                              "doc_id" : doc_id})
            output_list = []
            for attachment in list_of_attachments:
                attachment["base_url"] = schema.get("base_url", "")
                output = ATTACHMENT_TEMPLATE.format(**attachment)
                output_list.append(output)
            return ", \n".join(output_list)
                 
    return EMPTY_FILE_DESCRIPTION


def prepare_output(document, schema, export_type, slashed_json_pointer):
    
    input_data =  resolve_pointer(document,slashed_json_pointer, default='')

    if isinstance(input_data, basestring):
        return unicode(input_data)
    elif isinstance(input_data, list):
        return ", ".join([unicode(datum) for dataum in input_data])
    elif isinstance(input_data, dict):
        doc_id = document.get("id", "")        
        return prepare_file_output(input_data, schema, export_type, doc_id)
    return unicode(input_data)




def get_column( schema,documents, export_type):
    """Pull out the JSON pointer's value from the JSON doc"""

    slashed_json_pointer = "/%s" % schema["data"].replace(".", "/")
    return [prepare_output(document, schema, export_type, slashed_json_pointer) for document in documents]
           
            
        

def fix_column_types(df, schema):
    """Convert the columns in the dataframe to numerical types if appropriate"""
    for field in schema:
        if field.get("field_type", None) == "number":
            df[field["export_name"]] = df[field["export_name"]].astype(float)
        if field.get("field_type", None) == "integer":
            df[field["export_name"]] = df[field["export_name"]].astype(int)
        # if field.get("field_type", None) == "date":
        #     df[field["export_name"]] = pd.to_datetime(df[field["export_name"]], coerce=True)

def add_images_or_other_objects(col, column_schema, worksheet, objects, format, writer):
    """Generate and save images against a given column in the worksheet"""
    if column_schema.get("field_type", None) == "b64png":
        coldata = get_column(column_schema, objects, "xlsx")
        files = []
        for i, imagestring in enumerate(coldata):
            row_to_write = i+1
            worksheet.write_string(row_to_write, col, "")
            if len(imagestring) > 10:
                
                filename = '/tmp/' + str(objects[i]["id"]) +'.png'
                with open(filename, 'wb') as f:
                    f.write(base64.b64decode(imagestring))
                worksheet.set_row(row_to_write, 104, format)
                
                worksheet.insert_image(row_to_write, col, filename, {"x_offset": 6, "y_offset": 10, "x_scale": 1.29245283, "y_scale": 1.330188679})
                files.append(filename)
            else:
                worksheet.set_row(row_to_write, None, format)



class CBHCompoundBatchSerializer( Serializer):
    """Main serializer for CBHCompoundBatches, providing functions to output data as SDF or XLSX"""
    formats = ['json', 'jsonp', 'xml', 'yaml', 'html', 'csv', 'xlsx', 'sdf']
    content_types = {'json': 'application/json',
                     'jsonp': 'text/javascript',
                     'xml': 'application/xml',
                     'yaml': 'text/yaml',
                     'html': 'text/html',
                     'csv': 'text/csv',
                     'xlsx': 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                     'sdf': 'chemical/x-mdl-sdfile'}




    def to_xlsx(self, data, options=None):
        '''write excel file from the project by project data produced by the CBHCompoundbatchResource'''
        output = StringIO()

        project_records = self.to_simple(data, {})
        writer = pd.ExcelWriter('temp.xlsx', engine='xlsxwriter' , options={'encoding':'utf-8'})
        writer.book.filename = output
        if isinstance(data, dict):
            if data.get("traceback", False):
                # worksheet = workbook.add_worksheet("Error")
                # worksheet.write(0,0,data["traceback"])
                print data

        else:
            workbook = writer.book
            for projectsheet in project_records:
                # worksheet = workbook.add_worksheet(projectsheet["name"])
                # write_excel_headers(projectsheet["schema"], worksheet)
                jsondef = {}
                empty = False
                if len(projectsheet["objects"])== 0:
                    #Add a single default object to an empty dataframe to stop the bug in the writer
                    
                    empty = True
                else:
                    for col, column_schema in enumerate(projectsheet["schema"]):
                    #     write_excel_column(column_schema, projectsheet["objects"], worksheet, col)
                        jsondef[column_schema["export_name"]] = get_column(column_schema, projectsheet["objects"], "xlsx")
                    
                    df = pd.DataFrame(jsondef)
                    

                    df.fillna('', inplace=True)
                    df = df[[col["export_name"] for col in projectsheet["schema"]]]
                    if not empty:
                        fix_column_types(df, projectsheet["schema"])

                    df.to_excel(writer, sheet_name=projectsheet["name"], index=False,)
                    worksheet = writer.sheets[projectsheet["name"]]
                    
                    format = workbook.add_format()
                    format.set_text_wrap()
                    format.set_font_name('Cantarell')
                    
                    format2 = workbook.add_format({'bold': False})
                    format2.set_font_name('Cantarell')
                    format2.set_align('vcenter')
                    format2.set_align('center')
                    format2.set_bottom(2)
                    worksheet.set_row(0, 30, format2)
                    for col, schem in enumerate(projectsheet["schema"]):
                        worksheet.write(0,col, schem["export_name"] )
                    for col, column_schema in enumerate(projectsheet["schema"]):
                        add_images_or_other_objects(col, 
                                                        column_schema, 
                                                        worksheet,
                                                        projectsheet["objects"], 
                                                        format,
                                                        writer)    
                    widths = coreparser.get_widths(df)
                    for index, width in enumerate(widths):
                        if width > 150:
                            width = 150
                        elif width < 15:
                            width = 15
                        worksheet.set_column(index, index, width, format)
            writer.save()
        
       
        # for fi in list(set(files)):
        #     os.remove(fi)
        return output.getvalue()



    def to_sdf(self, data, options=None):
        ''' SDF converter'''
        project_records = self.to_simple(data, {})

        if isinstance(data, dict):
            if data.get("traceback", False):
                # worksheet = workbook.add_worksheet("Error")
                # worksheet.write(0,0,data["traceback"])
                print data

        mol_strings = []

        for projectsheet in project_records:
            # worksheet = workbook.add_worksheet(projectsheet["name"])
            # write_excel_headers(projectsheet["schema"], worksheet)
            
            jsondef = {}
            empty = False
            if len(projectsheet["objects"])== 0:
                #Add a single default object to an empty dataframe to stop the bug in the writer
                
                empty = True
            else:
                sdf_cols = [column_schema for column_schema in projectsheet["schema"] if column_schema["data"] != "image"]
                for col, column_schema in enumerate(sdf_cols):
                    jsondef[column_schema["export_name"]] = get_column(column_schema, projectsheet["objects"], "sdf")
                
                df = pd.DataFrame(jsondef)
                df.fillna('', inplace=True)
                df = df[[col["export_name"] for col in sdf_cols]]

                row_iterator = df.iterrows()
                headers = list(df)
                for index, row in row_iterator:
                    ctab = projectsheet["objects"][index]["ctab"]
                    if not ctab:
                        ctab = "\n\n\nM  END\n"
                    mol_str  = ctab.replace("RDKit          2D\n", "Generated by ChemBio Hub ChemiReg http://chembiohub.ox.ac.uk/chemireg\n").split("END")[0]
                    properties = []
                    for field in headers:
                        try:
                            properties.append(
                                SDF_TEMPLATE.format(
                                    **{"name": unicode(field), "value": unicode(row[field])}
                                )
                            )
                        except Exception as e:
                            import traceback
                            import sys
                            print(traceback.format_exception(*sys.exc_info()))


                    mol_strings.append("".join([mol_str, "END\n"] + properties))
        return "$$$$\n".join(mol_strings) + "$$$$"








