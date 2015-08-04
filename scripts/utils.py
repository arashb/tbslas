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
import parser
import pp

################################################################################
# GLOBALS
################################################################################
HOSTNAME      = socket.gethostname()
TIMESTR       = time.strftime("%Y%m%d-%H%M%S")
RESULT_TAG_HEADER  = '#TBSLAS-HEADER: '
RESULT_TAG         = '#TBSLAS-RESULT: '

################################################################################
# COMMANDLINE ARGUMENTS
################################################################################
USAGE = 'USAGE: ./python PROGRAM <mpi-num-processes> <omp-num-threads> <num-steps=5>'
print sys.argv
if len(sys.argv) < 3:
    print USAGE
    sys.exit()
if len(sys.argv) >= 3:
    MPI_TOTAL_NUM_PORCESSES = int(sys.argv[1])
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
OUTPUT_PREFIX   = SCRIPT_ID+'-np'+str(MPI_TOTAL_NUM_PORCESSES).zfill(5)+'-'+TIMESTR

def prepare_environment(output_prefix):
    global TBSLAS_RESULT_DIR_PREFIX
    RESULT_DIR = PWD
    if 'stampede' in HOSTNAME:
        RESULT_DIR = os.environ['SCRATCH']
    elif 'maverick' in HOSTNAME:
        RESULT_DIR = os.environ['WORK']
    print "STORING OUTPUT IN: " + RESULT_DIR
    TBSLAS_RESULT_DIR_PREFIX = os.path.join(RESULT_DIR, output_prefix)

def determine_command_prefix(mpi_num_procs):
    if 'stampede' in HOSTNAME:
        return ['ibrun', '-np', str(mpi_num_procs), 'tacc_affinity']
    elif 'maverick' in HOSTNAME:
        return ['ibrun', '-np', str(mpi_num_procs), 'tacc_affinity']
    else:
        return ['mpirun', '-np', str(mpi_num_procs)]

def analyse_command_output(output, \
                           file_dt, file_pr, file_out, file_pp, \
                           PRINT_RSLT_HEADER, PRINT_PRFL_HEADER):
    for line in output:
        file_out.write(line)
        # CATCH RESULTS HEADER
        if line.startswith(RESULT_TAG_HEADER) and PRINT_RSLT_HEADER:
            li = line.replace(RESULT_TAG_HEADER, '')
            file_dt.write(li)
            sys.stdout.write(li)
        # CATCH RESULTS DATA
        if line.startswith(RESULT_TAG):
            li = line.replace(RESULT_TAG, '')
            file_dt.write(li)
            sys.stdout.write(li)

    # PARSE PROFILE OUTPUT
    mydoc = parser.pdoc(output)
    mydoc.print_me(file_pr)
    pp.post_process_profile_data(mydoc, file_pp, PRINT_PRFL_HEADER);

def execute_commands(cmds, id):
    id = SCRIPT_ID+'-'+id
    sys.stdout.write("##############################\n")
    sys.stdout.write("# "+id+"\n")
    sys.stdout.write("##############################\n")
    PRINT_RSLT_HEADER = True
    PRINT_PRFL_HEADER = True

    # open output files
    if not os.path.exists(TBSLAS_RESULT_DIR_PREFIX):
        os.makedirs(TBSLAS_RESULT_DIR_PREFIX)
    file_dt = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.data'), 'w')
    file_pr = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.profile'), 'w')
    file_pp = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.profile.pp'), 'w')

    # output current git revision
    revision = subprocess.check_output(["git", "describe"])
    file_dt.write('# REVISION: ' + revision)
    file_pr.write('# REVISION: ' + revision)
    file_pp.write('# REVISION: ' + revision)

    # execute generated commands
    for counter, cmd in cmds.iteritems():
        out_dir_name = \
          os.path.join(TBSLAS_RESULT_DIR_PREFIX,'{0}-cmd{1:03}'.format(id,counter))
        out_file_name = out_dir_name+'.out'
        file_out = open(out_file_name, 'w')
        os.environ['TBSLAS_RESULT_DIR'] = out_dir_name
        os.makedirs(out_dir_name)

        # output command
        cmd_msg = '# CMD: ' +  ' '.join(cmd) + '\n'

        sys.stdout.write('# ============================================================================================================================================\n')
        sys.stdout.write(cmd_msg)

        file_out.write('# ============================================================================================================================================\n')
        file_out.write(cmd_msg)

        file_dt.write('# ============================================================================================================================================\n')
        file_dt.write(cmd_msg)

        file_pr.write('# ============================================================================================================================================\n')
        file_pr.write(cmd_msg)

        file_pp.write('# ============================================================================================================================================\n')
        file_pp.write (cmd_msg)

        # execute command
        p = subprocess.Popen(cmd,                    \
                             shell=False,            \
                             stdout=subprocess.PIPE, \
                             stderr=subprocess.STDOUT)

        # analyse command
        analyse_command_output(p.stdout.readlines(),                \
                               file_dt, file_pr, file_out, file_pp, \
                               PRINT_RSLT_HEADER, PRINT_PRFL_HEADER)
        PRINT_RSLT_HEADER = False
        PRINT_PRFL_HEADER = False

        # flush output
        file_dt.flush()
        file_pr.flush()
        file_pp.flush()
        file_out.close()
    file_dt.close()
    file_pr.close()
    file_pp.close()
