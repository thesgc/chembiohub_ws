# -*- coding: utf-8 -*-

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




def write_excel_headers(schema, worksheet):
    for column_no, field in enumerate(schema):
        worksheet.write(0, column_no, field["knownBy"])
    


def get_column( schema,documents):
    slashed_json_pointer = "/%s" % schema["data"].replace(".", "/")
    

    return [unicode(resolve_pointer(document,slashed_json_pointer, default=''))
            for  document in documents]
        

def fix_column_types(df, schema):
    for field in schema:
        if field.get("field_type", None) == "number":
            df[field["knownBy"]] = df[field["knownBy"]].astype(float)
        if field.get("field_type", None) == "integer":
            df[field["knownBy"]] = df[field["knownBy"]].astype(int)
        if field.get("field_type", None) == "date":
            df[field["knownBy"]] = pd.to_datetime(df[field["knownBy"]], coerce=True)

def add_images_or_other_objects(col, column_schema, worksheet, objects, format, writer):
    if column_schema.get("field_type", None) == "b64png":
        coldata = get_column(column_schema, objects)
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
        # for fi in list(set(files)):
        #     os.remove(fi)                      
            


class CBHCompoundBatchSerializer( Serializer):
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
        '''write excel file here'''
        output = StringIO()

        project_records = self.to_simple(data, {})
        writer = pd.ExcelWriter('temp.xlsx', engine='xlsxwriter' ,datetime_format='dd mmmm yyyy')
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
                        jsondef[column_schema["knownBy"]] = get_column(column_schema, projectsheet["objects"])
                    
                    df = pd.DataFrame(jsondef)
                    

                    df.fillna('', inplace=True)
                    df = df[[col["knownBy"].encode("ascii", "ignore") for col in projectsheet["schema"]]]
                    if not empty:
                        fix_column_types(df, projectsheet["schema"])

                    df.to_excel(writer, sheet_name=projectsheet["name"], index=False,)
                    worksheet = writer.sheets[projectsheet["name"]]
                    
                    format = workbook.add_format()
                    format.set_text_wrap()
                    format.set_align('vcenter')
                    format2 = workbook.add_format({'bold': False})

                    format2.set_align('vcenter')
                    format2.set_align('center')
                    format2.set_bottom(2)
                    worksheet.set_row(0, 30, None)
                    for col, schem in enumerate(projectsheet["schema"]):
                        worksheet.write(0,col, schem["knownBy"], format2)
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
        
            #workbook.close()
        

        # 

        # writer = pd.ExcelWriter('temp.xlsx', engine='xlsxwriter')
        # writer.book.filename = output
        # df.to_excel(writer, sheet_name='Sheet1', index=False)
        # workbook = writer.book
        # format = workbook.add_format()
        # worksheet = writer.sheets['Sheet1']
       
                



        # # make the UOx ID and SMILES columns bigger
        # # BUG - can't set column format until pandas 0.16
        # # https://github.com/pydata/pandas/issues/9167
        # for index, width in enumerate(widths):
        #     if width > 150:
        #         width = 150
        #     elif width < 15:
        #         width = 15
        #     worksheet.set_column(index, index, width)
        # writer.save()
        # for fi in list(set(files)):
        #     os.remove(fi)
        return output.getvalue()



    def to_sdf(self, data, options=None):
        '''Convert to SDF'''
        data = self.to_simple(data, {})
        try:
            if data.get("traceback", False):
                raise ImmediateHttpResponse(BadRequest(json.dumps(data.data)))
        except AttributeError:
            pass
        mols = []
        index = 0
        options = options or {}
        try:
            exp_json = json.loads(data.get('export', []))
        except:
            raise ImmediateHttpResponse(BadRequest(data))

        df = pd.DataFrame(exp_json)
        df.fillna('', inplace=True)
        cols = df.columns.tolist()
        # now for the list we have in the order we have it, move the columns by name
        # this way you end up with your core fields at the start and custom
        # fields at the end.
        ordered_fields = ['UOx ID', 'Project']
        headers = data.get('headers', {})
        ordered_fields += headers["custom_fields"]
        ordered_fields += headers["uncurated_fields"]
        for idx, item in enumerate(ordered_fields):
            cols.insert(idx, cols.pop(cols.index(item)))

        # reindex the dataframe
        df = df.ix[:, cols]
        # pull data back out of dataframe to put into rdkit tools

        row_iterator = df.iterrows()
        headers = list(df)
        mol_strings = []
        for index, row in row_iterator:
            # Simple string based SDF formatter
            mol_str = row['ctab'].replace(
                "RDKit          2D\n", "Generated by ChemBio Hub ChemiReg http://chembiohub.ox.ac.uk/chemireg\n").split("END")[0]
            properties = []
            for field in headers:
                if field != "ctab":
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








