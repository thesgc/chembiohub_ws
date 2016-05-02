from behave import given, when, then
#from django.contrib.auth.models import User



@when(u'I refresh the user object')
def step(context):
    context.user = User.objects.get(pk=context.user.pk)


