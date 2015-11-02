#*************************************************************************
#Copyright (C) 2015 by Arash Bakhtiari
#You may not use this file except in compliance with the License.
#You obtain a copy of the License in the LICENSE file.

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.
#*************************************************************************

################################################################################
# SYSTEM IMPORT
################################################################################
import subprocess
import socket
import time
import sys
import os

################################################################################
# GLOBALS
################################################################################
TIMESTR    = time.strftime("%Y%m%d-%H%M%S")
HOSTNAME   = socket.gethostname()
CMD_FFTW3  = ['module', 'load', 'fftw3']
CMD_PYTHON = ['module', 'load', 'python']

def submit_job(job_id, num_nodes, num_proces, num_threads, total_time):
    cmd_list = [];
    if 'stampede' in HOSTNAME:
        CMD_JOB =['sbatch', \
                  '-N'+str(num_nodes),\
                  '-n'+str(num_proces), \
                  '-p', 'normal', \
                  '-o '+job_id + TIMESTR + '.out', \
                  '--time='+str(total_time), \
                  '-J', job_id, \
                  './.run_python.sh', job_id, str(num_proces), str(num_threads)]
        cmd_list.extend([CMD_JOB])
    elif 'maverick' in HOSTNAME:
        CMD_JOB =['sbatch', \
                  '-N'+str(num_nodes),\
                  '-n'+str(num_proces), \
                  '-p', 'vis', \
                  '-o '+job_id + TIMESTR + '.out', \
                  '--time='+str(total_time), \
                  '-J', job_id, \
                  './.run_python.sh', job_id, str(num_proces), str(num_threads)]
        cmd_list.extend([CMD_JOB])
    elif 'zico' in HOSTNAME:
        CMD_JOB = ['./.run_python.sh', job_id, str(num_proces), str(num_threads)]
        cmd_list.extend([CMD_FFTW3, CMD_JOB])
    else:
        CMD_JOB = ['./.run_python.sh', job_id, str(num_proces), str(num_threads)]
        cmd_list.extend([CMD_JOB])
    # execute commands
    # print cmd_list
    for cmd in cmd_list:
        print cmd
        p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            sys.stdout.write(line)
