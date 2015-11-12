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

################################################################################
# GLOBALS
################################################################################
TIMESTR       = time.strftime("%Y%m%d-%H%M%S")
RESULT_TAG_HEADER  = '#TBSLAS-HEADER: '
RESULT_TAG         = '#TBSLAS-RESULT: '

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
TBSLAS_EXAMPLES_DIR      = os.path.join(TBSLAS_DIR, "examples/")
TBSLAS_EXAMPLES_BIN_DIR  = os.path.join(TBSLAS_EXAMPLES_DIR, "bin/")

################################################################################
# EXECUTION COMMAND
################################################################################
SCRIPT_ID       = sys.argv[0].replace('.py', '').replace('./','')

def parse_args():
    USAGE = 'USAGE: ./python PROGRAM <mpi-num-processes> <omp-num-threads>'
    print sys.argv
    if len(sys.argv) < 3:
        print USAGE
        sys.exit()
    if len(sys.argv) >= 3:
        mpi_num_procs = int(sys.argv[1])
        omp_num_threads = int(sys.argv[2])
    return (mpi_num_procs, omp_num_threads)

def get_output_prefix(num_proces):
    return SCRIPT_ID+'-np'+str(num_proces).zfill(5)+'-'+TIMESTR

def get_result_dir_prefix():
    hostname      = socket.gethostname()
    mpi_num_procs, omp_num_threads = parse_args()
    output_prefix = get_output_prefix(mpi_num_procs)
    tbslas_result_dir = os.environ['TBSLAS_RESULT_DIR']
    return (tbslas_result_dir, output_prefix)

def compile_code():
    print '--> compiling code ...'
    PWD = os.environ['PWD']
    os.chdir(TBSLAS_EXAMPLES_DIR)
    # execute command
    CMD_CLEAN = ['make', 'clean']
    CMD_MAKE = ['make', '-j']
    cmd_list = [CMD_CLEAN, CMD_MAKE]
    for cmd in cmd_list:
        try:
            proc = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            # print proc
        except subprocess.CalledProcessError as err:
            print err.output
            sys.exit()
    os.chdir(PWD)

def determine_command_prefix(mpi_num_procs, offset=0):
    hostname = socket.gethostname()
    if 'stampede' in hostname:
        return ['ibrun', '-n', str(mpi_num_procs), '-o', str(offset), 'tacc_affinity']
    elif 'maverick' in hostname:
        return ['ibrun', '-n', str(mpi_num_procs), '-o', str(offset), 'tacc_affinity']
    else:
        return ['mpirun', '-np', str(mpi_num_procs)]

def analyse_command_output(output, fout, fdata, fprof,
                           PRINT_RSLT_HEADER, PRINT_PRFL_HEADER,\
                           pp_func = None, fpp= None):
    print '--> analysing command output ...'
    for line in output:
        fout.write(line)
        # CATCH RESULTS HEADER
        if line.startswith(RESULT_TAG_HEADER) and PRINT_RSLT_HEADER:
            li = line.replace(RESULT_TAG_HEADER, '')
            fdata.write(li)
            sys.stdout.write(li)
        # CATCH RESULTS DATA
        if line.startswith(RESULT_TAG):
            li = line.replace(RESULT_TAG, '')
            fdata.write(li)
            sys.stdout.write(li)
    # PARSE PROFILE OUTPUT
    mydoc = parser.pdoc(output)
    mydoc.print_me(fprof)
    # POST PROCESSING
    if pp_func and fpp:
        pp_func(mydoc, fpp, PRINT_PRFL_HEADER);
        # pp.pp_profile_data(mydoc, fpp, PRINT_PRFL_HEADER);

def execute_commands(cmds, id, pp_func = None):
    id = SCRIPT_ID+'-'+id
    sys.stdout.write("##############################\n")
    sys.stdout.write("# "+id+"\n")
    sys.stdout.write("##############################\n")
    PRINT_RSLT_HEADER = True
    PRINT_PRFL_HEADER = True
    # open output files
    tbslas_result_dir, output_prefix = get_result_dir_prefix()
    TBSLAS_RESULT_DIR_PREFIX = os.path.join(tbslas_result_dir, output_prefix)
    if not os.path.exists(TBSLAS_RESULT_DIR_PREFIX):
        os.makedirs(TBSLAS_RESULT_DIR_PREFIX)
    flist = []
    # reporter output
    fdata = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.data'), 'w')
    flist.append(fdata)
    # profiling output
    fprof = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.prof'), 'w')
    flist.append(fprof)
    # post processing output
    fpp = None
    if pp_func:
        fpp = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.prof.pp'), 'w')
        flist.append(fpp)
    # output env metadata (git revision, host, ...)
    for f in flist:
        f.write('# REVISION: ' + subprocess.check_output(["git", "describe"]))
        f.write('# HOST: ' + socket.gethostname()+'\n')
    # execute generated commands
    for counter, cmd in cmds.iteritems():
        out_dir_name = \
          os.path.join(TBSLAS_RESULT_DIR_PREFIX,'{0}-cmd{1:03}'.format(id,counter))
        out_file_name = out_dir_name+'.out'
        fout = open(out_file_name, 'w')
        os.environ['TBSLAS_RESULT_DIR'] = out_dir_name
        os.makedirs(out_dir_name)
        # output command
        cmd_msg = '# CMD '+str(counter)+' : ' +  ' '.join(cmd) + '\n'
        sys.stdout.write('# ------------------------------\n')
        sys.stdout.write("--> storing output in: " + out_dir_name +" \n")
        sys.stdout.write("--> executing command ... \n")
        sys.stdout.write(cmd_msg)
        fout.write('# ------------------------------\n')
        fout.write(cmd_msg)
        for f in flist:
            f.write('# ------------------------------\n')
            f.write(cmd_msg)
        # execute command
        p = subprocess.Popen(cmd,                    \
                             shell=False,            \
                             stdout=subprocess.PIPE, \
                             stderr=subprocess.STDOUT)
        # analyse command
        analyse_command_output(p.stdout.readlines(),\
                               fout, fdata, fprof,\
                               PRINT_RSLT_HEADER, PRINT_PRFL_HEADER,\
                               pp_func, fpp)
        PRINT_RSLT_HEADER = False
        PRINT_PRFL_HEADER = False
        # flush output
        [f.flush() for f in flist]
        fout.close()
    [f.flush() for f in flist]
    # reset the env. variable
    os.environ['TBSLAS_RESULT_DIR'] = tbslas_result_dir
