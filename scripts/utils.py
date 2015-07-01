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
RESULT_TAG_HEADER  = '#TBSLAS-HEADER: '
RESULT_TAG_LIST    = ['#TBSLAS-RESULT: ']
PROFILE_TAG_HEADER = 't_min'
PROFILE_TAG_LIST   = ['+-RunSemilag', '+-RunFMM', '+-EvalTree', \
                      '+-MortonId', '+-ScatterIndex', '+-ScatterForward', '+-Evaluation', '+-ScatterReverse', 'TRG_CNT']
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
    elif 'maverick' in HOSTNAME:
        RESULT_DIR = os.environ['WORK']
    print "STORING OUTPUT IN: " + RESULT_DIR
    TBSLAS_RESULT_DIR_PREFIX = os.path.join(RESULT_DIR, output_prefix)

def determine_command_prefix():
    if 'stampede' in HOSTNAME:
        return ['ibrun', 'tacc_affinity']
    elif 'maverick' in HOSTNAME:
        return ['ibrun', 'tacc_affinity']
    else:
        return ['mpirun', '-n', str(MPI_NUM_PROCESS)]

def analyse_command_output(output, \
                           file_dat, file_prf, file_out, \
                           PRINT_RSLT_HEADER, PRINT_PRFL_HEADER):
    for line in output:
        file_out.write(line)
        # CATCH RESULTS HEADER
        if line.startswith(RESULT_TAG_HEADER) and PRINT_RSLT_HEADER:
            li = line.replace(RESULT_TAG_HEADER, '')
            file_dat.write(li)
            sys.stdout.write(li)
        # CATCH RESULTS DATA
        for result_tag in RESULT_TAG_LIST:
            if line.startswith(result_tag):
                li = line.replace(result_tag, '')
                file_dat.write(li)
                sys.stdout.write(li)
        # CATCH PROFILE HEADER
        if PROFILE_TAG_HEADER in line and PRINT_PRFL_HEADER:
            li = 'MODULE' + line[6:]
            sys.stdout.write(li)
            file_prf.write('# ============================================================================================================================================\n')
            file_prf.write(li)
            file_prf.write('# ============================================================================================================================================\n')
        # CATCH PROFILE DATA
        for profile_tag in PROFILE_TAG_LIST:
            if profile_tag in line:
                li = line.replace("+-", '')
                sys.stdout.write(li)
                file_prf.write(li)
                # file_prf.write('# ----------------------------------------------------------------------\n')

def execute_commands(cmd_args, id):
    id = SCRIPT_ID+'-'+id
    sys.stdout.write("##############################\n")
    sys.stdout.write("# "+id+"\n")
    sys.stdout.write("##############################\n")
    PRINT_RSLT_HEADER = True
    PRINT_PRFL_HEADER = True
    # metadata and convergence results
    if not os.path.exists(TBSLAS_RESULT_DIR_PREFIX):
        os.makedirs(TBSLAS_RESULT_DIR_PREFIX)
    file_dat = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.dat'), 'w')
    file_prf = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.prf'), 'w')
    revision = subprocess.check_output(["git", "describe"])
    file_dat.write('# REVISION: '+ revision)
    file_prf.write('# REVISION: '+ revision)
    for counter, args in cmd_args.iteritems():
        out_dir_name = \
          os.path.join(TBSLAS_RESULT_DIR_PREFIX,'{0}-cmd{1:03}'.format(id,counter))
        out_file_name = out_dir_name+'.out'
        file_out = open(out_file_name, 'w')
        command_message = "COMMAND: " +  str(args) + '\n'
        sys.stdout.write(command_message)
        os.environ['TBSLAS_RESULT_DIR'] = out_dir_name
        os.makedirs(out_dir_name)
        # execute command
        p = subprocess.Popen(determine_command_prefix() + args, \
                             shell=False,                       \
                             stdout=subprocess.PIPE,            \
                             stderr=subprocess.STDOUT)
        analyse_command_output(p.stdout.readlines(), \
                               file_dat, file_prf, file_out, \
                               PRINT_RSLT_HEADER, PRINT_PRFL_HEADER)
        PRINT_RSLT_HEADER = False
        PRINT_PRFL_HEADER = False
        file_dat.flush()
        file_prf.flush()
        file_out.close()
    file_dat.close()
    file_prf.close()
