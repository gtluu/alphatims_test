# views.py
from flask import abort, jsonify, render_template, request, redirect, url_for, make_response, send_from_directory
import uuid

from app import app
import os
import glob
import json
import tarfile

UPLOAD_FOLDER = "./data"

@app.route('/heartbeat', methods=['GET'])
def heartbeat():
    return "{}"

@app.route('/convert', methods=['POST'])
def convert():
    # Generate UUID.
    job_uuid = str(uuid.uuid4().hex)

    # Creating outputing directory
    temp_dir = os.path.join(UPLOAD_FOLDER, job_uuid)
    os.makedirs(temp_dir) # Making sure the folder exists

    # Get tarball and save to server.
    uploaded_data = request.files['data']
    uploaded_data_path = os.path.join(UPLOAD_FOLDER, job_uuid, 'upload.tar.gz')
    uploaded_data.save(uploaded_data_path)

    # Decompress uploaded data.
    with tarfile.open(uploaded_data_path) as tarball:
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tarball, temp_dir)

    # Build TIMSCONVERT command.
    run_script = "/app/timsconvert/bin/run.py"
    input_file = glob.glob(os.path.join(temp_dir, 'data', '*.d'))[0]
    output_file = 'output'
    # hard code exclude mobility for now
    cmd = 'python {} --input {} --outfile {} --exclude_mobility'.format(run_script, input_file, output_file)

    # Run TIMSCONVERT
    os.system(cmd)

    # Tar output files.
    output_tar = os.path.join(temp_dir, 'output.tar.gz')
    if os.path.exists(os.path.join(temp_dir, 'data', 'output.mzML')):
        with tarfile.open(output_tar, 'w:gz') as newtar:
            newtar.add(os.path.join(temp_dir, 'data', 'output.mzML'), 'output.mzML')
    elif os.path.exists(os.path.join(temp_dir, 'data', 'output.imzML')) and os.path.exists(os.path.join(temp_dir, 'data', 'output.ibd')):
        with tarfile.open(output_tar, 'w:gz') as newtar:
            newtar.add(os.path.join(temp_dir, 'data', 'output.imzML'), 'output.imzML')
            newtar.add(os.path.join(temp_dir, 'data', 'output.ibd'), 'output.ibd')

    # Send files to client.
    return send_from_directory(temp_dir, 'output.tar.gz')
