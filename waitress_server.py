import logging

from waitress import serve
import Executor
print('Server started at port 8080')
logger = logging.getLogger('waitress')
logger.setLevel(logging.INFO)
logger.info("Hellow world")
Executor.app.secret_key='12345'
serve(Executor.app, host='0.0.0.0', port=8080, threads=100)
