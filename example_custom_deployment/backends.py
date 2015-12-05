from django.conf import settings
from django.contrib.auth.models import User, check_password, Group


def test_login_on_external_system(username, password):
    """
    replace return True with the appropriate login code for your system
    """
    return True



class CustomAuthBackend(object):
    """
    Authenticates a user via an external system
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
                g = Group.objects.get_or_create(name="Standard Permissions Group")
                g.user_set.add(user)
                g.save()
            return user
        return None

    def get_user(self, user_id):
        try:
            return User.objects.get(pk=user_id)
        except User.DoesNotExist:
            return None