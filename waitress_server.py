import logging

from waitress import serve
import Executor
print('Server started at port 8080')
Executor.app.logger.setLevel(logging.DEBUG)
serve(Executor.app, host='0.0.0.0', port=8080, threads=8)
