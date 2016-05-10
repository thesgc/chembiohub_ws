from behave import given, when, then

@when(u'I refresh the user object')
def step(context):
    from cbh_core_model.models import User
    context.user = User.objects.get(pk=context.user.pk)
