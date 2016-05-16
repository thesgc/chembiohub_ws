import os
import datetime
from settings import FLOWJS_PATH, FLOWJS_EXPIRATION_DAYS


def chunk_upload_to(instance, filename):
    """
    Save chunk to the right path and filename based in is number
    """
    return os.path.join(FLOWJS_PATH, instance.filename)


