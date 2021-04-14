import os
import getpass
import subprocess
from platform import platform
import sys
import argparse
import pickle
import json # For dealing with metadata
from datafed.CommandLib import API

parser = argparse.ArgumentParser(description='Create datasets for DataFed')
parser.add_argument('--id', dest='id', action='store',
                   help='Collection id to create data within')
parser.add_argument('--f', dest='file', action='store',
                   help='Filename')
parser.add_argument('--n', dest='name', action='store',
                  help='Name of dataset')
parser.add_argument('--m', dest='metadata', action='store',
                  help='Metadata file')
args = parser.parse_args()
filename = args.file
file_title=args.name
global_coll_id = args.id
df_api = API()
pkl_file = args.metadata



with open(pkl_file, 'rb') as f:
    metadata = pickle.load(f)
rec_msg = df_api.dataCreate(title = file_title,
                            alias = '',
                            metadata=json.dumps(metadata),
                            parent_id=global_coll_id,
                                )
rec_id = rec_msg[0].data[0].id
#Use as pathname the path and name of the file you wish to move from CADES to DataFed
pput_msg = df_api.dataPut(rec_id, filename, wait=False)
