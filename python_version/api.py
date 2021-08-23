import os
import time
import traceback
from datetime import datetime

from flask import Flask
from flask import request, jsonify

from Common.Service.benchmark import Benchmark

api = Flask(__name__)
b = Benchmark()
rs = ResourceService()

api.wsgi_app = Middleware(api.wsgi_app)
# Reset stats for each request
api.before_request(f=b.reset_time)
# Clean temporal data after every request
api.after_request(f=rs.clean_temporal_images)

