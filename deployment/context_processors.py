from django.conf import settings

def static_paths(request):
    print settings.ENV_NAME
    # if settings.ENV_NAME  != "dev":
    #     return {"site_prefix" : "/%s/" % settings.ENV_NAME}
    # else:
    return {"site_prefix" : "/../" }