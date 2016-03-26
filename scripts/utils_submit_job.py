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
import uuid

################################################################################
# LOCAL IMPORT
################################################################################
import utils

################################################################################
# GLOBALS
################################################################################
HOSTNAME   = socket.getfqdn()


def get_supermuc_batch_file (nnodes, nprocs, nthreads, myqueue, job_id, TIMESTR, total_time):
        print os.environ['HOME']

#template job script
        commands=\
"""#! /usr/bin/ksh
#@ shell = /usr/bin/ksh
#@ job_type = MPICH
#@ initialdir={TBSLAS_DIR}/scripts
#@ job_name = {JOB_ID}
#@ class = {QUEUE}
#@ node_usage = not_shared
#@ wall_clock_limit = {TOTAL_TIME}
#@ network.MPI = sn_all,,us,,
#@ notification = never
#@ output = {JOB_ID}.$(jobid).out
#@ error =  {JOB_ID}.$(jobid).err

#@ node = {NUM_NODES}
#@ total_tasks = {NUM_PROCS}
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

module load python/2.7_anaconda
module load git
module load fftw/serial/3.3
module load gcc/4.9

export TBSLAS_RESULT_DIR={TBSLAS_RESULT_DIR}
export TBSLAS_DIR={TBSLAS_DIR}

module unload mpi.ibm
module load mpi.intel
./.run_python.sh {JOB_id} {NUM_PROCS} {NUM_THREADS}\n"""

        output_dir='./'
        # if not os.path.exists(output_dir):
        #         os.makedirs(output_dir)

        file_name= os.path.join(output_dir, job_id[:-3]+'_'+TIMESTR+'.cmd')
        file_handler = open(file_name,"w")

        commands=commands.format(JOB_ID=job_id[:-3],
                                 JOB_id=job_id,
                                 QUEUE=myqueue,
                                 TOTAL_TIME=total_time,
                                 NUM_NODES=nnodes,
                                 NUM_PROCS=nprocs,
                                 NUM_THREADS=nthreads,
                                 TBSLAS_RESULT_DIR=os.environ['TBSLAS_RESULT_DIR'],
                                 TBSLAS_DIR=os.environ['TBSLAS_DIR'])

        file_handler.write(commands)
        file_handler.close()

        return file_name

def submit_job(job_id, num_nodes, num_procs, num_threads, total_time=None, queue=None):
    print '--> submit job ' + job_id + ' ...'
    TIMESTR    = time.strftime("%Y%m%d-%H%M%S-")+str(uuid.uuid4())

    cmd_list = [];
    if not total_time:
        total_time = '01:00:00'
    if 'stampede' in HOSTNAME:          # TACC Stampede cluster
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
    elif 'maverick' in HOSTNAME:        # TACC Maverick cluster
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
    elif 'cos.lrz.de' in HOSTNAME:      # LRZ linux cluster
        if not queue:
            queue = 'mpp2'
        CMD_JOB =['sbatch', \
                  '-N'+str(num_nodes),\
                  '-n'+str(num_procs), \
                  '-M', queue, \
                  '-o'+job_id+'_'+TIMESTR + '.out', \
                  '--time='+str(total_time), \
                  '-J', job_id, \
                  './.run_python.sh', job_id, str(num_procs), str(num_threads)]
        cmd_list.extend([CMD_JOB])
    elif 'sm.lrz.de' in HOSTNAME:      # LRZ SuperMUC cluster
        if not queue:
            queue = 'micro'

        jobfile = get_supermuc_batch_file(num_nodes, num_procs, num_threads, queue, job_id, TIMESTR, total_time)
        CMD_JOB = ['llsubmit',str(jobfile)]
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
    if len(sys.argv) < 5:
        print USAGE
        sys.exit()
    if len(sys.argv) >= 5:
        job_id      = sys.argv[1]
        num_nodes   = int(sys.argv[2])
        num_procs   = int(sys.argv[3])
        num_threads = int(sys.argv[4])
    total_time=None
    if len(sys.argv) >= 6:
        total_time  = sys.argv[5]
    queue = None
    if len(sys.argv) >= 7:
        queue = sys.argv[6]
    submit_job(job_id, num_nodes, num_procs, num_threads, total_time, queue)
