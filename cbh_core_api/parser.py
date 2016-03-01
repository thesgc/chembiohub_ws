
from pandas import ExcelFile
import xlrd
from copy import copy

def get_widths(df):
    '''Returns a list of the maximum string widths of a the columns in a dataframe in characters'''
    widths = []
    for col in df.columns.tolist():
        col = unicode(col)
        titlewidth = len(col)
        try:
            w = df[col].values.astype(unicode).str.len().max()
            if w > titlewidth:
                widths.append(int(w*1.2))
            else:
                widths.append(int(titlewidth * 1.2))
        except:
            widths.append(int(titlewidth * 1.2))
    return widths


def is_true(item):
    if str(item).lower() in ["y", "true", "yes"]:
        return True
    else:
        return False


def get_custom_field_config(filename, sheetname):
    '''Early example of the import of a custom field config based upon a list of field names'''
    xls = ExcelFile(filename)
    data = xls.parse(sheetname, index_col=None, na_values=[''])
    data.columns = ["name", "required", "description"]
    data["required"] = data["required"].apply(is_true)

    data = data.fillna('')
    mydata = [{key: unicode(value) for key, value in point.items()} for point in data.T.to_dict().values()]
    return mydata


def get_key_from_field_name(name):
    '''Returns the field names used by elasticsearch'''
    return unicode(name).replace(u" ", u"__space__")


def get_sheetnames(filename):
    '''Returns a list of the sheet names in an Excel file'''
    xls = xlrd.open_workbook(filename, on_demand=True)
    return xls.sheet_names()


def get_sheet(filename, sheetname):
    '''Extracts a list of dicts from a worksheet of an Excel file along with the
    column names, data types and maximum widths'''
    xls = ExcelFile(filename)
    data = xls.parse(sheetname, index_col=None, na_values=[''])
    data = data.fillna('')
    orig_cols = tuple(data.columns)
    replace = [get_key_from_field_name(column) for column in data.columns]
    data.columns = replace
    types = copy(data.dtypes)
    for col in replace:
        data[col] = data[col].values.astype(unicode)
    return (data.T.to_dict().values(), orig_cols, types, get_widths(data))


