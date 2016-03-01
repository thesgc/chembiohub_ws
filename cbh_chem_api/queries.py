from tastypie.resources import Resource
from tastypie import fields
from django.forms import Form

SORT_DIRECTION_CHOICES = (("asc", "Ascending"),
                           ("desc", "Decending"),
                           ("none", "No Sort") )

QUERY_TYPE_CHOICES = (("any", "Any in this list"),
                        ("exactly", "Equals"),
                        ("phrase", "Contains phrase starting with"),
                        ("starts_with", "Starts with"),
                        ("ends_with", "Ends with"),
                        ("between", "Between"),
                        ("greater_than", "Greater than"),
                        ("less_than", "Greater than"),)

class CBHQueryForm(Form):
    name = forms.CharField(max_length=100)
    query_type = forms.ChoiceField(max_length=20, choices=QUERY_TYPE_CHOICES)
    equals = forms.CharField(max_length=100)
    maximum = forms.CharField(max_length=20)
    mimimum = forms.CharField(max_length=20)
    sort_direction = forms.ChoiceField(max_length=20, choices=SORT_DIRECTION_CHOICES)




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
