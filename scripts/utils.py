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

import subprocess
import socket
import time
import sys
import os

################################################################################
# GLOBALS
################################################################################
HOSTNAME      = socket.gethostname()
TIMESTR       = time.strftime("%Y%m%d-%H%M%S")
HEADER_TAG    = '#TBSLAS-HEADER: '
RESULT_TAG    = '#TBSLAS-RESULT: '
################################################################################
# COMMANDLINE ARGUMENTS
################################################################################
USAGE = 'USAGE: ./python PROGRAM <mpi-num-processes> <omp-num-threads> <num-steps=5>'
print sys.argv
if len(sys.argv) < 3:
    print USAGE
    sys.exit()
if len(sys.argv) >= 3:
    MPI_NUM_PROCESS = int(sys.argv[1])
    OMP_NUM_THREADS = int(sys.argv[2])
################################################################################
# ENVIRONMENT VARIABLES
################################################################################
try:
    TBSLAS_DIR = os.environ['TBSLAS_DIR']
except KeyError as e:
    print "Environment variable {0} is not set.".format(e)
    sys.exit()
################################################################################
# DIRECTORIES
################################################################################
PWD = os.environ['PWD']
TBSLAS_EXAMPLES_DIR = os.path.join(TBSLAS_DIR, "examples/")
TBSLAS_EXAMPLES_BIN_DIR = os.path.join(TBSLAS_EXAMPLES_DIR, "bin/")
TBSLAS_RESULT_DIR_PREFIX = ''
################################################################################
# EXECUTION COMMAND
################################################################################
SCRIPT_ID       = sys.argv[0].replace('.py', '').replace('./','')
OUTPUT_PREFIX = SCRIPT_ID+'-'+TIMESTR
def prepare_environment(output_prefix):
    global TBSLAS_RESULT_DIR_PREFIX
    RESULT_DIR = PWD
    if 'stampede' in HOSTNAME:
        RESULT_DIR = os.environ['SCRATCH']
    print "STORING OUTPUT IN: " + RESULT_DIR
    TBSLAS_RESULT_DIR_PREFIX = os.path.join(RESULT_DIR, output_prefix)

def determine_command_prefix():
    global HOSTNAME
    global MPI_NUM_PROCESS
    if 'stampede' in HOSTNAME:
        return ['ibrun', 'tacc_affinity']
    else:
        return ['mpirun', '-n', str(MPI_NUM_PROCESS)]

def analyse_command_output(output, file_dat, file_out, PRINT_HEADER):
    for line in output:
        file_out.write(line)
        if line.startswith(HEADER_TAG) and PRINT_HEADER:
            li = line.replace(HEADER_TAG, '')
            file_dat.write(li)
            sys.stdout.write(li)
        elif line.startswith(RESULT_TAG):
            li = line.replace(RESULT_TAG, '')
            file_dat.write(li)
            sys.stdout.write(li)

def execute_commands(cmd_args, id):
    id = SCRIPT_ID+'-'+id
    sys.stdout.write("##############################\n")
    sys.stdout.write("# "+id+"\n")
    sys.stdout.write("##############################\n")
    PRINT_HEADER = True
    # metadata and convergence results
    if not os.path.exists(TBSLAS_RESULT_DIR_PREFIX):
        os.makedirs(TBSLAS_RESULT_DIR_PREFIX)
    file_dat = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.dat'), 'w')
    revision = subprocess.check_output(["git", "describe"])
    file_dat.write('# REVISION: '+ revision)
    for counter, args in cmd_args.iteritems():
        out_file_name = os.path.join(TBSLAS_RESULT_DIR_PREFIX,'{i}-cmd{x}.out'.format(i=id,x=counter))
        out_dir_name = os.path.join(TBSLAS_RESULT_DIR_PREFIX,'{i}-cmd{x}'.format(i=id,x=counter))
        file_out = open(out_file_name, 'w')
        command_message = "COMMAND: " +  str(args) + '\n'
        #sys.stdout.write(command_message)
        os.environ['TBSLAS_RESULT_DIR'] = out_dir_name
        os.makedirs(out_dir_name)
        # execute command
        p = subprocess.Popen(determine_command_prefix() + args, \
                             shell=False,                       \
                             stdout=subprocess.PIPE,            \
                             stderr=subprocess.STDOUT)
        analyse_command_output(p.stdout.readlines(), file_dat, file_out, PRINT_HEADER)
        PRINT_HEADER = False
        file_dat.flush()
        file_out.close()
    file_dat.close()
