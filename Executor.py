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
app = Flask(__name__)
import sys
import shutil
from flask import Response
import enhancerGeneDistanceBased
import eQTLAanalysis
import mergeMethods
from multiprocessing import Process
import datetime
import threading

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
        if (current_time - created_time).total_seconds() > 480 or (current_time - modified_time).total_seconds() > 480:
            if os.path.isdir(file_path):
                shutil.rmtree(file_path)
                logging.info("deleted folder ",file_path)
            elif os.path.isfile(file_path) or os.path.islink(file_path):
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

def executeFunc(a,organ):
    df1 = pd.read_csv('mapping.csv')
    organFiles = df1[df1['Organ'] == organ]
    chiapet = organFiles['peekachu']
    eqtl = organFiles['gtex']
    eqtlHelp = organFiles['getxHelper']
    temp = "temp/"
    chiaFile = temp+'tempChiapet'+str(uuid.uuid4())+'.bed'
    distanceFile = temp+'tempFinaldistance'+str(uuid.uuid4())+'.bed'
    eqtlFile = temp+'tempFinaleQTL' + str(uuid.uuid4()) + '.bed'

    # ray.get([chiaPetAnalysis.startPoint.remote(a,chiapet.values[0],chiaFile),
    #          enhancerGeneDistanceBased.startPoint.remote(a,distanceFile),
    #          eQTLAanalysis.startPoint.remote(a,eqtl.values[0],eqtlHelp.values[0],eqtlFile)])

    chiaPetAnalysis.startPoint(a,chiapet.values[0],chiaFile)
    enhancerGeneDistanceBased.startPoint(a,distanceFile)
    eQTLAanalysis.startPoint(a,eqtl.values[0],eqtlHelp.values[0],eqtlFile)
    imagesFileName = temp+'images-'+str(uuid.uuid4())
    os.mkdir(imagesFileName[0:len(imagesFileName)])
    mergeMethods.startPoint(chiaFile,distanceFile,eqtlFile,imagesFileName+'/')
    os.remove(chiaFile)
    os.remove(distanceFile)
    os.remove(eqtlFile)
    shutil.make_archive(imagesFileName, 'zip', imagesFileName)
    imagesFileName = imagesFileName[len(temp)-1:len(imagesFileName)]
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
                filename = secure_filename(file.filename)
                file.save(filename)
                images = executeFunc(filename, organ)
                return render_template('final.html', filename=images)
                # return sendFile(images)
            else:
                raise Exception('file name is not properly formatted')

    except Exception as e:
        abort(make_response(jsonify(message=str(e)), 400))



  #returns the username
  # return 'Username: %s' % username



# In[4]:





# In[ ]:


# os.system('python3 mergeMethods.py '+a)

