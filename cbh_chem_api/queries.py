from tastypie.resources import Resource
from tastypie import fields
from django.forms import Form




def get_single_nested_term_query( fieldname, term):
    pass











class CBHQueryField(Resource):
    name = fields.CharField(max_length=100)
    query_type = fields.CharField(max_length=20)
    equals = fields.CharField(max_length=100)
    maximum = fields.CharField(max_length=20)
    mimimum = fields.CharField(max_length=20)
    sort_direction = fields.CharField(max_length=20)

    class Meta:
        form = FormValidation(form_class=CBHQueryForm)


class CBHSearch(Resource):
    query_fields = fields.ToManyField(CBHQueryField, full=True)
