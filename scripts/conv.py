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
# LOCAL IMPORT
################################################################################
import submit_job as sj
import utils

################################################################################
# PARSE ARGUMENTS
################################################################################
USAGE = 'USAGE: ./python PROGRAM <num-nodes> <mpi-num-processes> <omp-num-threads> <job-duration>'
print sys.argv

num_nodes   = 1
num_procs   = 1
num_threads = 1
total_time  = '04:00:00'

if len(sys.argv) < 4:
    print USAGE
    sys.exit()
if len(sys.argv) >= 4:
    num_nodes   = int(sys.argv[1])
    num_procs   = int(sys.argv[2])
    num_threads = int(sys.argv[3])
if len(sys.argv) >= 5:
    total_time = sys.argv[4]

################################################################################
# COMPILE CODE
################################################################################
utils.compile_code()

################################################################################
# SUBMIT CONVERGENCE JOBS
################################################################################

job_list = [\
    # 'advection',\
    # 'diffusion',\
    # 'advdiff',\
    'zalesak',\
    # 'cubic',\
    ]
for job_id in job_list:
    sj.submit_job('conv-'+job_id+'.py', num_nodes, num_procs, num_threads, total_time)
