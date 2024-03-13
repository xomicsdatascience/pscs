from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

limiter = Limiter(get_remote_address, default_limits=["200 per day", "50 per hour"], storage_uri="mongodb://127.0.0.1:27017")