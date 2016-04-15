# -*- coding: utf-8 -*-
"""
Serializer functions required to export a list of fields for a project
"""
import cStringIO

from tastypie.serializers import Serializer

import pandas as pd
import numpy as np

from tastypie.exceptions import ImmediateHttpResponse, BadRequest
import json
from django.conf import settings

def get_field_name_from_key(key):
    """probably deprecated"""
    return key.replace(u"__space__", u" ")


def get_key_from_field_name(name):
    """probably deprecated"""
    return name.replace(u" ", u"__space__")


class CustomFieldXLSSerializer(Serializer):

    ''' Code for preparing an Excel summary of the custom fields for the given project '''
    formats = ['json', 'jsonp', 'xlsx']
    content_types = {'json': 'application/json',
                     'jsonp': 'text/javascript',
                     'xlsx': 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'}




    def to_xlsx(self, data, options=None):
        """Translate a project's fields into xlsx"""

        try:
            #Annoying way in which errors during serialization then cause tastypie to
            #try to serialize the error response as XLSX which obviously doesn't work so
            #We catch any response with a traceback in it and ensure it gets send straight to the front end
            if data.data.get("traceback", False):
                raise ImmediateHttpResponse(BadRequest(json.dumps(data.data)))
        except AttributeError:
            pass

        output = cStringIO.StringIO()
        try:
            exp_json = data.get('custom_field_config', None)
        except AttributeError:
            exp_json = None

        if exp_json is None:

            exp_json = data.data["project_data_fields"]
            exp_json = [field.data for field in exp_json]

        cleaned_data = []

        # need to manipulate the dataset which is used to apply to dataframe
        # date fields do not have allowed values but do have specified data
        # ranges
        for field in exp_json:

            # is it a date field? add the date ranges to the allowed values
            # column
            if(field['field_type'] == 'date'):
                field['field_type'] = 'Date'
                field['allowed_values'] = 'Valid Date'
            elif(field['field_type'] == 'uiselecttags'):
                field['field_type'] = 'Multiple select'
                field[
                    'placeholder'] = 'Select one or more of the Allowed Values, separated by a comma'
            elif(field['field_type'] == 'percentage'):
                field['field_type'] = 'Percentage'
                field['allowed_values'] = '0-100%'
            elif(field['field_type'] == 'text'):
                field['field_type'] = 'Plain Text'
                field[
                    'allowed_values'] = 'Plain Text. There may be a character length restriction on this field.'

        df = pd.DataFrame(
            exp_json, columns=['name', 'field_type', 'description', 'allowed_values'])
        # human readable titles
        df.rename(columns={'name': 'Name', 'field_type': 'Data Type', 'description':
                           'Description', 'allowed_values': 'Allowed Values'}, inplace=True)

        # deal with empty fields
        df.fillna('', inplace=True)

        # autosize column widths setup
        widths = []
        for col in df.columns.tolist():
            col = str(col)
            titlewidth = len(col)
            try:
                w = df[col].astype(unicode).str.len().max()
                if w > titlewidth:
                    widths.append(int(w*1.2))
                else:
                    widths.append(int(titlewidth * 1.2))
            except:
                widths.append(int(titlewidth * 1.2))

        writer = pd.ExcelWriter('temp.xlsx', engine='xlsxwriter')
        writer.book.filename = output
        df2 = pd.DataFrame(data=np.zeros((0, len(exp_json))), columns=[
                           field["name"] for field in exp_json])
        df2.to_excel(writer, sheet_name='Sheet1', index=False)

        df.to_excel(writer, sheet_name='Sheet2', index=False)

        workbook = writer.book
        format = workbook.add_format()
        worksheet = writer.sheets['Sheet2']
        format.set_text_wrap()
        # make the UOx ID and SMILES columns bigger
        # BUG - can't set column format until pandas 0.16
        # https://github.com/pydata/pandas/issues/9167
        for index, width in enumerate(widths):
            if width > 150:
                width = 150
            elif width < 15:
                width = 15
            worksheet.set_column(index, index, width)

        worksheet = writer.sheets['Sheet1']
        for index, field in enumerate(exp_json):
            width = len(field["name"])*1.5
            worksheet.set_column(index, index, width)
        writer.save()
        return output.getvalue()






