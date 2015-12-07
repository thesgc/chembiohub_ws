from django.conf import settings
from django.contrib.auth.models import User, check_password, Group
from django.contrib.auth.backends import ModelBackend

def test_login_on_external_system(username, password):
    """
    replace return True with the appropriate login code for your system
    """
    return True



class CustomAuthBackend(ModelBackend):
    """
    Authenticates a user via an external system
    We subclass the model backend to ensure the get user and all the permissions functions exist
    """

    def authenticate(self, username=None, password=None):
        login_valid = test_login_on_external_system(username, password)
        if login_valid:
            try:
                user = User.objects.get(username=username)
            except User.DoesNotExist:
                # Create a new user. Note that we can set password
                # to anything, because it won't be checked; the password
                # from the external system will
                user = User(username=username)
                
                if username in settings.ADMIN_USERS:
                    user.is_staff = True
                    user.is_superuser = True

                user.save()
                g, created = Group.objects.get_or_create(name="Standard Permissions Group")
                user.groups.add(g)
                user.save()
            return user
        return None
