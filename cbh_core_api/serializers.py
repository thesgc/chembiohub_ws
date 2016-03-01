# -*- coding: utf-8 -*-

import cStringIO

from tastypie.serializers import Serializer

import pandas as pd
import numpy as np

from tastypie.exceptions import ImmediateHttpResponse, BadRequest
import json

def get_field_name_from_key(key):
    return key.replace(u"__space__", u" ")


def get_key_from_field_name(name):
    return name.replace(u" ", u"__space__")


class CustomFieldXLSSerializer(Serializer):

    ''' COde for preparing an Excel summary of the custom fields for the given project '''
    formats = ['json', 'jsonp', 'xlsx']
    content_types = {'json': 'application/json',
                     'jsonp': 'text/javascript',
                     'xlsx': 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'}

    def to_xlsx(self, data, options=None):
        try:
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


class CustomFieldsSerializer(CustomFieldXLSSerializer):
    pass


class ResultsExportXLSSerializer(Serializer):

    ''' COde for preparing an Excel summary of the results from a particular search query '''
    formats = ['json', 'jsonp', 'xlsx']
    content_types = {'json': 'application/json',
                     'jsonp': 'text/javascript',
                     'xlsx': 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'}

    def to_xlsx(self, data, options=None):

        output = cStringIO.StringIO()

        exp_json = self.to_simple(data, {})

        if exp_json.get("traceback", False):
            raise ImmediateHttpResponse(json.dumps(exp_json))
        cleaned_data = []

        chembl_data = []

        col_titles = {}

        # do your stuff here
        # maybe just get all the l3 data out? group by their l2 uri, so the
        # columns might match up?
        for result, val in exp_json.iteritems():

            if(result == 'hits'):
                for k, v in val.iteritems():
                    if (k == 'hits'):
                        for item in v:

                            src = item['_source']['l3']['project_data']
                            for sk, sv in src.iteritems():
                                readable_sk = get_field_name_from_key(sk)
                                col_titles[sk] = readable_sk
                            srcl2 = item['_source']['l2']['project_data']
                            srcl1 = item['_source']['l1']['project_data']
                            srcl0 = item['_source']['l0']['project_data']
                            # print(srcl2['Title'])

                            src['l2'] = srcl2.get('Title', None)
                            src['l1'] = srcl1.get('Title', None)
                            src['l0'] = srcl0.get('Title', None)

                            col_titles['l2'] = "Assay"
                            col_titles['l1'] = "Study"
                            col_titles['l0'] = "Project"

                            cleaned_data.append(src)

                            # now chembl data!
                            if(item['_source']['l3'].get('chembl', None)):
                                csrc = item['_source']['l3']['chembl']
                                # are the chembl fields stored in that format?
                                # for sk, sv in csrc.iteritems():
                                #readable_sk = get_field_name_from_key(sk)
                                #col_titles[sk] = readable_sk

                                chembl_data.append(csrc)

        # for k, v in l3.iteritems():
        #    k = get_field_name_from_key(k)

        #df = pd.DataFrame(exp_json, columns=['name', 'field_type', 'description', 'allowed_values'])
        df = pd.DataFrame(cleaned_data)

        # human readable titles
        # now doing human readable titles via get_field_name_from_key
        df.rename(columns=col_titles, inplace=True)

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
        #df2 = pd.DataFrame(data=np.zeros((0,len(exp_json))), columns=[field["name"] for field in exp_json])
        #df2.to_excel(writer, sheet_name='Sheet1', index=False)

        df.to_excel(writer, sheet_name='Sheet1', index=False)

        if chembl_data:
            cdf = pd.DataFrame(chembl_data)
            cdf.fillna('', inplace=True)
            cdf.to_excel(writer, sheet_name='Sheet2', index=False)

        workbook = writer.book
        format = workbook.add_format()
        # worksheet = writer.sheets['Sheet2']

        worksheet = writer.sheets['Sheet1']

        # auto adjust column widths
        for index, width in enumerate(widths):
            if width > 150:
                width = 150
            elif width < 15:
                width = 15
            worksheet.set_column(index, index, width)
        writer.save()
        return output.getvalue()
