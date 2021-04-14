import os
import getpass
import subprocess
from platform import platform
import sys

#Verify Datafed installation
try:
    # This package is not part of anaconda and may need to be installed.
    from datafed.CommandLib import API
except ImportError:
    print('datafed not found. Installing from pip.')
    subprocess.call([sys.executable, "-m", "pip", "install", 'datafed'])
    from datafed.CommandLib import API

from datafed import version as df_ver

if not df_ver.startswith('1.1.0'):
    print('Attempting to update DataFed.')
    subprocess.call([sys.executable, "-m", "pip", "install", '--upgrade', 'datafed'])
    print('Please restart the python kernel or upgrade manually to V 1.1.0:1 if you are repeatedly seeing this message via'
          '\n\tpip install --upgrade datafed')
else:
    df_api = API()
    #print('Success! You have DataFed: ' + df_ver)
# Verify user authentication
if not df_api.getAuthUser():
    print('You have not authenticated into DataFed Client')


# Check default Globus endpoint
if not df_api.endpointDefaultGet():
    endpoint = 'cades#CADES-OR'
    df_api.endpointDefaultSet(endpoint)



#print('Your default Globus Endpoint in DataFed is:\n' + df_api.endpointDefaultGet())
# Test the endpoint
dget_resp = df_api.dataGet('d/35437908',
                           '/~/',
                           wait=True)
#print(dget_resp)
if dget_resp[0].task[0].status == 3:
    #os.remove('35437908.md5sum')
    sys.exit(0)
else:
    if dget_resp[0].task[0].msg == "globus connect offline":
        print('You need to activate your Globus Endpoint and/or ensure Globus Connect Personal is running.\n'
              'Please visit https://globus.org to activate your Endpoint')
        sys.exit(1)
    elif dget_resp[0].task[0].msg == "permission denied":
        print('Globus does not have write access to this directory. \n'
              'If you are using Globus Connect Personal, ensure that this notebook runs within'
              'one of the directories where Globus has write access. You may consider moving this'
              'notebook to a valid directory or add this directory to the Globus Connect Personal settings')
        sys.exit(1)
    else:
        NotImplementedError('Get in touch with us or consider looking online to find a solution to this problem:\n' + dget_resp[0].task[0].msg)
        sys.exit(1)
