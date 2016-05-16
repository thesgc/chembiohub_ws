from django.contrib import admin
from cbh_chembl_id_generator.models import CBHPlugin


class CBHPluginAdmin(admin.ModelAdmin):
    list_display = ('full_function_name', 'plugin_type', 'input_json_path')

admin.site.register(CBHPlugin, CBHPluginAdmin)
