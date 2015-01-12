#*************************************************************************
#Copyright (C) 2014 by Arash Bakhtiari
#You may not use this file except in compliance with the License.
#You obtain a copy of the License in the LICENSE file.

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.
#*************************************************************************
import os
import subprocess
import math
import time
import socket
import sys

from collections import OrderedDict
################################################################################
# GLOBALS
################################################################################
if len(sys.argv) < 5:
    raise
MPI_NUM_PROCESS = int(sys.argv[1])
OMP_NUM_THREADS = int(sys.argv[2])
TOL_NUM_DIGITS_INIT = int(sys.argv[3])
TOL_NUM_DIGITS_FINAL = int(sys.argv[4])
TOL_NUM_STEPS = TOL_NUM_DIGITS_FINAL - TOL_NUM_DIGITS_INIT + 1 
HOSTNAME = socket.gethostname()
TIMESTR  = time.strftime("%Y%m%d-%H%M%S")
OUTPUT_PREFIX = 'conv-'+TIMESTR

################################################################################
# ENVIRONMENT VARIABLES
################################################################################
try:
    TBSLAS_DIR = os.environ['TBSLAS_DIR']
except KeyError as e:
    print "Environment variable {0} is not set.".format(e)
    sys.exit()

# os.environ['OMP_NUM_THREADS'] = str(OMP_NUM_THREADS)
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
PROGRAM  = "zalesak"
EXEC     = os.path.join(TBSLAS_EXAMPLES_BIN_DIR, PROGRAM)

def prepare_environment():
    global TBSLAS_RESULT_DIR_PREFIX
    RESULT_DIR = PWD
    if 'stampede' in HOSTNAME:
        RESULT_DIR = os.environ['WORK']
    TBSLAS_RESULT_DIR_PREFIX = os.path.join(RESULT_DIR, OUTPUT_PREFIX)

def determine_command_prefix():
    if 'stampede' in HOSTNAME:
        return ['ibrun', 'tacc_affinity']
    else:
        return ['mpirun', '-n', str(MPI_NUM_PROCESS)]

def generate_commands():
    # generate a dictionary data type of commands
    commands = OrderedDict()
    dt = 0.01 #math.pi/30
    tn = 1
    tol = math.pow(0.1,float(TOL_NUM_DIGITS_INIT))
    # tol_list = [math.pow(0.1,x) for x in range(TOL_NUM_DIGITS_INIT,TOL_NUM_DIGITS_FINAL+1)]
    for counter in range(0,TOL_NUM_STEPS):
    # for tol in tol_list:
        ARGS    = ['-N', '8', '-tol', str(tol), '-dt', str(dt), '-tn', str(tn), '-omp', str(OMP_NUM_THREADS)]
        cmd = determine_command_prefix() + [EXEC] + ARGS
        # save command
        commands[tol] = cmd
        dt  = dt*0.5
        tol = tol*0.5
        tn  = 2*tn
    return commands

def execute_commands(commands):
    output = []
    for tolerance, command in commands.iteritems():
        command_message = "COMMAND: " +  str(command) + '\n'
        sys.stdout.write(command_message)
        os.environ['TBSLAS_RESULT_DIR'] = os.path.join(TBSLAS_RESULT_DIR_PREFIX,'tol_{x:.2e}'.format(x=tolerance))
        os.makedirs(os.environ['TBSLAS_RESULT_DIR'])
        # execute command
        p = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            if line.startswith('TOL:'):
                sys.stdout.write(line)
            output.append(line)
    return output

def analyse_output(output):
    # metadata and convergence results
    res_file_meta = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX,OUTPUT_PREFIX+'.metadata'), 'w')
    res_file_conv = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX,OUTPUT_PREFIX+'.out'), 'w')
    for line in output:
        res_file_meta.write(line)
        if line.startswith('TOL:'):
            res_file_conv.write(line)
    res_file_conv.close()
    res_file_meta.close()

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    commands = generate_commands()
    prepare_environment()
    output   = execute_commands(commands)
    analyse_output(output)
