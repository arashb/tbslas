#!/usr/bin/python

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

################################################################################
# GLOBALS
################################################################################
MPI_NUM_PROCESS = sys.argv[1]
OMP_NUM_THREADS = sys.argv[2]

HOSTNAME = socket.gethostname()
TIMESTR  = time.strftime("%Y%m%d-%H%M%S")

################################################################################
# ENVIRONMENT VARIABLES
################################################################################
try:
    TBSLAS_DIR = os.environ['TBSLAS_DIR']
except KeyError as e:
    print "Environment variable {0} is not set.".format(e)
    sys.exit()

os.environ['OMP_NUM_THREADS'] = OMP_NUM_THREADS
################################################################################
# DIRECTORIES
################################################################################
TBSLAS_EXAMPLES_DIR = os.path.join(TBSLAS_DIR, "examples/")
TBSLAS_EXAMPLES_BIN_DIR = os.path.join(TBSLAS_EXAMPLES_DIR, "bin/")

################################################################################
# EXECUTION COMMAND
################################################################################
PROGRAM  = "zalesak"
EXEC     = os.path.join(TBSLAS_EXAMPLES_BIN_DIR, PROGRAM)

def determine_command_prefix():
    if 'stampede' in HOSTNAME:
        return ['ibrun', 'tacc_affinity']
    elif 'zico' in HOSTNAME:
        return ['mpirun', '-n', str(MPI_NUM_PROCESS)]

    return []

def generate_commands():
    commands = []
    dt = 0.1047 #math.pi/30
    tn = 1
    toli = 3
    tolf = 4
    tol_list = [math.pow(0.1,x) for x in range(toli,tolf+1)]
    for tol in tol_list:
        ARGS    = ['-N', '8', '-tol', str(tol), '-dt', str(dt), '-tn', str(tn)]
        cmd = determine_command_prefix() + [EXEC] + ARGS
        # save command
        commands.append(cmd)
        dt = dt*0.5
        tn = 2*tn
    return commands

def execute_commands(commands):
    output = []
    for command in commands:
        command_message = "COMMAND: " +  str(command) + '\n'
        print(command_message)
        # output = output +  [command_message]
        # execute command
        p = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            print line
            output.append(line)
    return output

def analyse_output(output):
    f = open('conv-'+TIMESTR+'.out','w')
    for line in output:
        if line.startswith('TOL:'):
            f.write(line)
    f.close()

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    commands = generate_commands()
    output   = execute_commands(commands)
    analyse_output(output)
