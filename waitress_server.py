from waitress import serve
import Executor
print('Server started at port 8080')
serve(Executor.app, host='0.0.0.0', port=8080)
