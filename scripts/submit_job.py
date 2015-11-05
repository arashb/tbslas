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

def submit_job(job_id, num_nodes, num_procs, num_threads, total_time, queue=None):
    cmd_list = [];
    if 'stampede' in HOSTNAME:
        if not queue:
            queue = 'normal'
        CMD_JOB =['sbatch', \
                  '-N'+str(num_nodes),\
                  '-n'+str(num_procs), \
                  '-p', queue, \
                  '-o'+ job_id + TIMESTR + '.out', \
                  '--time='+str(total_time), \
                  '-J', job_id, \
                  './.run_python.sh', job_id, str(num_procs), str(num_threads)]
        cmd_list.extend([CMD_JOB])
    elif 'maverick' in HOSTNAME:
        if not queue:
            queue = 'vis'
        CMD_JOB =['sbatch', \
                  '-N'+str(num_nodes),\
                  '-n'+str(num_procs), \
                  '-p', queue, \
                  '-o'+job_id+'_'+TIMESTR + '.out', \
                  '--time='+str(total_time), \
                  '-J', job_id, \
                  './.run_python.sh', job_id, str(num_procs), str(num_threads)]
        cmd_list.extend([CMD_JOB])
    elif 'zico' in HOSTNAME:
        CMD_JOB = ['./.run_python.sh', job_id, str(num_procs), str(num_threads)]
        cmd_list.extend([CMD_JOB])
    else:
        CMD_JOB = ['./.run_python.sh', job_id, str(num_procs), str(num_threads)]
        cmd_list.extend([CMD_JOB])
    # execute commands
    # print cmd_list
    for cmd in cmd_list:
        print cmd
        p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            sys.stdout.write(line)

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    USAGE = 'USAGE: python submit_job.py <job-id> <num-nodes> <mpi-num-processes> <omp-num-threads> <total-time>'
    print sys.argv
    if len(sys.argv) < 6:
        print USAGE
        sys.exit()
    if len(sys.argv) >= 6:
        job_id      = sys.argv[1]
        num_nodes   = int(sys.argv[2])
        num_procs   = int(sys.argv[3])
        num_threads = int(sys.argv[4])
        total_time  = sys.argv[5]
    queue = None
    if len(sys.argv) >= 7:
        queue = sys.argv[6]
    # print len(sys.argv)
    # print sys.argv
    submit_job(job_id, num_nodes, num_procs, num_threads, total_time, queue)
