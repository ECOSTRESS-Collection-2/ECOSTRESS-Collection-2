import base64
import json
import os
from pathlib import Path
from typing import Optional

import earthaccess


class CMRServerUnreachable(Exception):
    pass


_AUTH = None
def login() -> earthaccess.Auth:
    # Only login to earthaccess once
    global _AUTH
    if _AUTH is not None:
        return _AUTH

    if not "EARTHDATA_USERNAME" in os.environ or not "EARTHDATA_PASSWORD" in os.environ:
        raise CMRServerUnreachable("Missing environment variable 'EARTHDATA_USERNAME' or 'EARTHDATA_PASSWORD'")

    try:
        _AUTH = earthaccess.login(strategy="environment")
        return _AUTH
    except Exception as e:
        raise CMRServerUnreachable(e)
