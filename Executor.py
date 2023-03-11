#!/usr/bin/env python
# coding: utf-8

# In[4]:
import os
import uuid
import pandas as pd
from flask import Flask, flash, request, abort, make_response, jsonify
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

def executeFunc(a,organ):
    df1 = pd.read_csv('mapping.csv')
    organFiles = df1[df1['Organ'] == organ]
    chiapet = organFiles['peekachu']
    eqtl = organFiles['gtex']
    eqtlHelp = organFiles['getxHelper']

    chiaFile = 'tempChiapet'+str(uuid.uuid4())+'.bed'
    distanceFile = 'tempFinaldistance'+str(uuid.uuid4())+'.bed'
    eqtlFile = 'tempFinaleQTL' + str(uuid.uuid4()) + '.bed'
    chiaPetAnalysis.startPoint(a,chiapet.values[0],chiaFile)
    enhancerGeneDistanceBased.startPoint(a,distanceFile)
    eQTLAanalysis.startPoint(a,eqtl.values[0],eqtlHelp.values[0],eqtlFile)
    imagesFileName = 'images-'+str(uuid.uuid4())
    os.mkdir(imagesFileName[0:len(imagesFileName)])
    mergeMethods.startPoint(chiaFile,distanceFile,eqtlFile,imagesFileName+'/')
    os.remove(chiaFile)
    os.remove(distanceFile)
    os.remove(eqtlFile)
    shutil.make_archive(imagesFileName, 'zip', imagesFileName)
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
def sendFile(imagesFileName):
    path=''
    zipname = imagesFileName+'.zip'
    with open(os.path.join(path, zipname), 'rb') as f:
        data = f.readlines()
    os.remove(os.path.join(path, zipname))
    shutil.rmtree(os.path.join(path,imagesFileName))
    return Response(data, headers={
        'Content-Type': 'application/zip',
        'Content-Disposition': 'attachment; filename=%s;' % zipname
    })



def allowed_file(filename):
    ALLOWED_EXTENSIONS = ['bed', 'gz']
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/execute/<organ>',methods=['GET', 'POST'])
def show_user(organ):
    try:
        if request.method == 'POST':
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
                images = executeFunc(filename,organ)
                return sendFile(images)
            else:
                raise Exception('file name is not properly formatted')

    except Exception as e:
        abort(make_response(jsonify(message=str(e)), 400))




  #returns the username
  # return 'Username: %s' % username



# In[4]:





# In[ ]:


# os.system('python3 mergeMethods.py '+a)

