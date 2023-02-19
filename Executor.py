#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os
import pandas as pd
from flask import Flask, flash, request, abort, make_response, jsonify
from werkzeug.utils import secure_filename

app = Flask(__name__)
import sys
import shutil
from flask import Response

# In[ ]:


def executeFunc(a,organ):
    df1 = pd.read_csv('mapping.csv')
    organFiles = df1[df1['Organ'] == organ]
    chiapet = organFiles['peekachu']
    eqtl = organFiles['gtex']
    eqtlHelp = organFiles['getxHelper']
    os.system('python3 chiaPet-Analysis.py '+a +' '+ chiapet.values[0])
    os.system('python3 enhancer-gene-distanceBased.py '+a+' & python3 eQTL-Aanalysis.py '+a+' '+eqtl.values[0] +' '+eqtlHelp.values[0])
    os.system('python3 mergeMethods.py '+a)
    shutil.make_archive('file-output', 'zip', 'images')
    


# In[3]:
def sendFile():
    path=''
    zipname = 'file-output.zip'
    with open(os.path.join(path, zipname), 'rb') as f:
        data = f.readlines()
    os.remove(os.path.join(path, zipname))
    return Response(data, headers={
        'Content-Type': 'application/zip',
        'Content-Disposition': 'attachment; filename=%s;' % zipname
    })

ALLOWED_EXTENSIONS = ['bed','bed.gz']
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS
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
                executeFunc(filename,organ)
                return sendFile()
            else:
                raise Exception('file name is not properly formatted')

    except Exception as e:
        abort(make_response(jsonify(message=str(e)), 400))


    
   
  #returns the username
  # return 'Username: %s' % username



# In[4]:





# In[ ]:


# os.system('python3 mergeMethods.py '+a)

