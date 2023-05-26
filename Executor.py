#!/usr/bin/env python
# coding: utf-8
import logging
# In[4]:
import os
import uuid
import pandas as pd
from flask import Flask, flash, request, abort, make_response, jsonify, render_template, send_from_directory
from werkzeug.utils import secure_filename
import chiaPetAnalysis
# from Process import Process
from Process import Process

app = Flask(__name__)
app.logger.setLevel(logging.INFO)
import sys
import shutil
from flask import Response
import enhancerGeneDistanceBased
import eQTLAanalysis
import mergeMethods
import datetime
import threading
from pybedtools.helpers import BEDToolsError, cleanup
import pandas as pd

logger = logging.getLogger('waitress')
logger.setLevel(logging.INFO)

def my_scheduled_job():
    logging.info("Scheduler running")
    directory_path = 'temp/'

# get the current time
    current_time = datetime.datetime.now()

    # loop through all the files in the directory
    for file in os.listdir(directory_path):

        # get the full file path
        file_path = os.path.join(directory_path, file)

        # get the file's creation time and modified time
        created_time = datetime.datetime.fromtimestamp(os.path.getctime(file_path))
        modified_time = datetime.datetime.fromtimestamp(os.path.getmtime(file_path))

        # check if the file was created or modified more than 8 minutes ago
        if (current_time - modified_time).total_seconds() > 900 or (current_time - modified_time).total_seconds() > 900:
            if os.path.isdir(file_path):
                shutil.rmtree(file_path)
                logging.info("deleted folder ",file_path)
            elif os.path.isfile(file_path) or os.path.islink(file_path):
                if(file_path !='temp/tempr'):
                    os.remove(file_path)
                    logging.info("deleted file ", file_path)


def my_method():
    # your code here
    my_scheduled_job()

def run_method():
    print("started")
    my_method()
    threading.Timer(60, run_method).start()

run_method()

def target():
    raise ValueError('Something went wrong...')

def startPosition(df):
    return 'EH'+str(uuid.uuid4())
def executeFunc(a,organ):
    df1 = pd.read_csv('mapping.csv')
    organFiles = df1[df1['Organ'] == organ]
    chiapet = organFiles['peekachu']
    eqtl = organFiles['gtex']
    eqtlHelp = organFiles['getxHelper']
    temp = "temp/"
    # if not os.path.exists(temp):
    #     os.makedirs(temp)
    myuuid = str(uuid.uuid4())
    imagesFileName = temp + 'images-' + myuuid
    os.mkdir(imagesFileName[0:len(imagesFileName)])
    chiaFile = imagesFileName+'/tempChiapet'+myuuid+'.bed'
    distanceFile = imagesFileName+'/tempFinaldistance'+myuuid+'.bed'
    eqtlFile = imagesFileName+'/tempFinaleQTL' + myuuid + '.bed'
    df22 = pd.read_csv( a, sep='\t', header=None)
    tempEQTL1 = 'temp/rawBed' + str(uuid.uuid4()) + '.bed'
    if(len(df22.columns) < 4):
        df22[3] = df22.apply(startPosition, axis=1)
    else:
        df22 = df22[[0, 1,2,3]]
    df22 = df22.sort_values(by=[0,1, 2], ascending=True)
    new_header = df22.iloc[0]
    df22 = df22[1:]
    df22.columns = new_header
    df22.to_csv(tempEQTL1, sep="\t", index=False)
    p1 = Process(target=chiaPetAnalysis.startPoint,args=[tempEQTL1,chiapet.values[0],chiaFile])
    p2 = Process(target=enhancerGeneDistanceBased.startPoint,args=[tempEQTL1,distanceFile])
    p3 = Process(target=eQTLAanalysis.startPoint,args=[tempEQTL1,eqtl.values[0],eqtlHelp.values[0],eqtlFile])
    # try:
    p1.start()
    p2.start()
    p3.start()
    p1.join()
    p2.join()
    p3.join()
    if p1.exception or p2.exception or p3.exception:
        print(p1.exception)
        if  "chromomsome sort ordering" in str(p1.exception) or "chromomsome sort ordering" in str(p2.exception) or "chromomsome sort ordering" in str(p3.exception):
            raise Exception("The file provided does not belong to the associated organ")
        else:
            raise Exception("Something went wrong with the file, please check your upload file")


    # except BEDToolsError as e:
    #     print("error:", e, file=sys.stderr)
    #     os._exit(1)
    #     raise Exception("Something is wrong")





    # result_chiaPetAnalysis = pool.map_async(chiaPetAnalysis.startPoint(a,chiapet.values[0],chiaFile))
    # result_enhancerGeneDistanceBased = pool.map_async(enhancerGeneDistanceBased.startPoint(a,distanceFile))
    # result_eQTLAanalysis = pool.map_async(eQTLAanalysis.startPoint(a,eqtl.values[0],eqtlHelp.values[0],eqtlFile))
    # print(result_eQTLAanalysis,result_chiaPetAnalysis,result_enhancerGeneDistanceBased)

    # ray.get([chiaPetAnalysis.startPoint.remote(a,chiapet.values[0],chiaFile),
    #          enhancerGeneDistanceBased.startPoint.remote(a,distanceFile),
    #          eQTLAanalysis.startPoint.remote(a,eqtl.values[0],eqtlHelp.values[0],eqtlFile)])

    # chiaPetAnalysis.startPoint(a,chiapet.values[0],chiaFile)
    # enhancerGeneDistanceBased.startPoint(a,distanceFile)
    # eQTLAanalysis.startPoint(a,eqtl.values[0],eqtlHelp.values[0],eqtlFile)

    mergeMethods.startPoint(chiaFile,distanceFile,eqtlFile,imagesFileName+'/')

    shutil.make_archive(imagesFileName, 'zip', imagesFileName)
    imagesFileName = imagesFileName[len(temp)-1:len(imagesFileName)]
    os.remove(chiaFile)
    os.remove(distanceFile)
    os.remove(eqtlFile)
    return imagesFileName

def runInParallel(*fns):
  proc = []
  for fn in fns:
    p = Process(target=fn)
    p.start()
    proc.append(p)
  for p in proc:
    p.join()

# In[3]:
@app.route('/database_download/<filename>')
def database_download(filename):
    path='temp'

    zipname = filename+'.zip'
    with open(os.path.join(path, zipname), 'rb') as f:
        data = f.readlines()
    os.remove(os.path.join(path, zipname))
    shutil.rmtree(os.path.join(path,filename))
    return Response(data, headers={
        'Content-Type': 'application/zip',
        'Content-Disposition': 'attachment; filename=%s;' % zipname
    })


def allowed_file(filename):
    ALLOWED_EXTENSIONS = ['bed', 'gz']
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# @app.route('/database_download/<filename>')
# def database_download(filename):
#     return send_from_directory('database_reports', filename)
# @app.route('/execute/<organ>',methods=['GET', 'POST'])
# def show_user(organ):
#     try:
#         if request.method == 'POST':
#             # check if the post request has the file part
#             if 'file' not in request.files:
#                 flash('No file part')
#                 raise Exception('No file part')
#             file = request.files['file']
#             if file.filename == '':
#                 flash('No selected file')
#                 raise Exception('No Selected file')
#             if file and allowed_file(file.filename):
#                 filename = secure_filename(file.filename)
#                 file.save(filename)
#                 images = executeFunc(filename,organ)
#                 return render_template('final.html', filename=images)
#                 # return sendFile(images)
#             else:
#                 raise Exception('file name is not properly formatted')
#
#     except Exception as e:
#         abort(make_response(jsonify(message=str(e)), 400))

@app.route('/', methods=['GET','POST'])
def hello():
    try:
        if request.method =='GET':
            return render_template('index1.html')
        elif request.method == 'POST':
            organ = request.form['organ']
            if(organ == 'Select an organ'):
                raise Exception('No Selected organ')
            print(organ)
            # check if the post request has the file part
            if 'file' not in request.files:
                flash('No file part')
                raise Exception('No file part')
            file = request.files['file']
            if file.filename == '':
                flash('No selected file')
                raise Exception('No Selected file')
            if file and allowed_file(file.filename):
                filename = "temp/"+secure_filename(file.filename)
                file.save(filename)
                images = executeFunc(filename, organ)
                return render_template('final.html', filename=images)
                # return sendFile(images)
            else:
                raise Exception('file name is not properly formatted')

    except Exception as e:
        logger.error(str(e))
        raise Exception('file name is not properly formatted')
        abort(make_response(jsonify(message=str(e)), 400))



  #returns the username
  # return 'Username: %s' % username



# In[4]:





# In[ ]:


# os.system('python3 mergeMethods.py '+a)

