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
from collections import OrderedDict

################################################################################
# LOCAL IMPORT
################################################################################
import utils_parser as parser
import utils_pp as pp

################################################################################
# GLOBALS
################################################################################
TIMESTR       = time.strftime("%Y%m%d-%H%M%S")+str(time.time())
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

def get_output_prefix():
    # return SCRIPT_ID+'-np'+str(num_proces).zfill(5)+'-'+TIMESTR
    return SCRIPT_ID+'-'+TIMESTR#+'-'+str(uuid.uuid4())


def get_result_dir_prefix():
    hostname      = socket.gethostname()
    # mpi_num_procs, omp_num_threads = parse_args()
    output_prefix = get_output_prefix()
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
#     hostname = socket.gethostname()
    hostname   = socket.getfqdn()

    if 'stampede' in hostname:
        return ['ibrun', '-n', str(mpi_num_procs), '-o', str(offset), 'tacc_affinity']
    elif 'maverick' in hostname:
        return ['ibrun', '-n', str(mpi_num_procs), '-o', str(offset), 'tacc_affinity']
    elif 'nid' in hostname:
        return ['ibrun', '-n', str(mpi_num_procs), '-o', str(offset), 'tacc_affinity']
    elif 'sm.lrz.de' in hostname:
        return ['mpiexec', '-n', str(mpi_num_procs)]
    elif 'cos.lrz.de' in hostname:
        return ['mpiexec', '-n', str(mpi_num_procs)]
    else:
        return ['mpirun', '-np', str(mpi_num_procs)]

def generate_commands(
    # executable name
    prog,
    # number of points
    pn_list = None,
    # tree tolerance
    tl_list = None,
    # maximux tree depth
    dp_list = None,
    # cheybishev degree
    cq_list = None,
    # cubic interpolation
    ci_list = None,
    # upsampling factor in cubic interp
    uf_list = None,
    # number of processes
    np_list = None,
    # number of threads
    nt_list = None,
    # temporal resolution
    dt_list = None,
    # number of time steps
    tn_list = None,
    # save VTK files
    vs_list = None,
    # load-balancing (tree merge) method
    mg_list = None,
    # test number
    tt_list = None,
    # diffusivity for advection-diffusion solver
    di_list = None,
    # alpha value for advection-diffusion solver
    ea_list = None):

    num_steps = len(pn_list)
    EXEC = os.path.join(TBSLAS_EXAMPLES_BIN_DIR, prog)

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS = []
        if pn_list:
            ARGS.extend(['-N'     , str(pn_list[counter])])
        if tl_list:
            ARGS.extend(['-tol'   , str(tl_list[counter])])
        if dp_list:
            ARGS.extend(['-d'     , str(dp_list[counter])])
        if cq_list:
            ARGS.extend(['-q'     , str(cq_list[counter])])
        if uf_list:
            ARGS.extend(['-cuf'   , str(uf_list[counter])])
        if tn_list:
            ARGS.extend(['-tn'    , str(tn_list[counter])])
        if dt_list:
            ARGS.extend(['-dt'    , str(dt_list[counter])])
        if mg_list:
            ARGS.extend(['-merge' , str(mg_list[counter])])
        if nt_list:
            ARGS.extend(['-omp'   , str(nt_list[counter])])
        if vs_list:
            ARGS.extend(['-vsr'   , str(vs_list[counter])])
        if tt_list:
            ARGS.extend(['-test'   , str(tt_list[counter])])
        if di_list:
            ARGS.extend(['-diff'   , str(di_list[counter])])
        if ea_list:
            ARGS.extend(['-ea'   , str(ea_list[counter])])
        if ci_list[counter]:
            ARGS = ARGS + ['-cubic', '1']
        cmd_args[cmd_id] = determine_command_prefix(np_list[counter]) + [EXEC] + ARGS
        print "CMD {0}: {1} ".format(cmd_id, cmd_args[cmd_id])
        cmd_id = cmd_id + 1
    return cmd_args

def analyse_command_output(output, fout, fdata, fprof,
                           PRINT_RSLT_HEADER, PRINT_PRFL_HEADER,\
                           pp_func = None, fpp= None):
    print '--> analysing command output ...'
    li_header = []
    li_values = []
    for line in output:
        fout.write(line)
        # CATCH RESULTS HEADER
        if line.startswith(RESULT_TAG_HEADER):
            li_header = line.replace(RESULT_TAG_HEADER, '').rstrip('\n')
            continue
        # CATCH RESULTS DATA
        if line.startswith(RESULT_TAG):
            li_values = line.replace(RESULT_TAG, '').rstrip('\n')
            continue

    # PARSE PROFILE OUTPUT
    mydoc = parser.pdoc(output)
    mydoc.print_me(fprof)

    # ADD TIMING VALUES TO REPORT
    tvals = pp.get_time(mydoc)
    for key, val in tvals.iteritems():
         li_header += "{:>10}".format(key)
         li_values += "{:>10.2f}".format(val)

    if PRINT_RSLT_HEADER:
        fdata.write(li_header + "\n")
    fdata.write(li_values + "\n")

def execute_commands(cmds, id, pp_func = None):
    sys.stdout.write("##############################\n")
    sys.stdout.write("# "+id+"\n")
    sys.stdout.write("##############################\n")
    PRINT_RSLT_HEADER = True
    PRINT_PRFL_HEADER = True
    # open output files
    tbslas_result_dir, output_prefix = get_result_dir_prefix()
    TBSLAS_RESULT_DIR_PREFIX = tbslas_result_dir
    # TBSLAS_RESULT_DIR_PREFIX = os.path.join(tbslas_result_dir, output_prefix)
    if not os.path.exists(TBSLAS_RESULT_DIR_PREFIX):
        os.makedirs(TBSLAS_RESULT_DIR_PREFIX)
    flist = []
    # reporter output
    fdata = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.data'), 'w')
    flist.append(fdata)
    # profiling output
    fprof = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.prof'), 'w')
    flist.append(fprof)

    # commands list
    fcmds = open(os.path.join(TBSLAS_RESULT_DIR_PREFIX, id+'.cmds'), 'w')
    flist.append(fcmds)

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
    for cmd_id, cmd in cmds.iteritems():
        out_dir_name = \
          os.path.join(TBSLAS_RESULT_DIR_PREFIX,'{0}-cmd{1:03}'.format(id,cmd_id))

        os.environ['TBSLAS_RESULT_DIR'] = out_dir_name
        os.makedirs(out_dir_name)

        sys.stdout.write('# ------------------------------\n')
        sys.stdout.write("--> storing output in: " + out_dir_name +" \n")
        sys.stdout.write("--> executing command {0} ... \n".format(cmd_id))

        # output file
        out_file_name = out_dir_name+'.out'
        fout = open(out_file_name, 'w')

        # output command
        cmd_msg = '# CMD '+str(cmd_id)+' : ' +  ' '.join(cmd) + '\n'
        fout.write('# ------------------------------\n')
        fout.write(cmd_msg)
        fcmds.write(cmd_msg)
        fcmds.flush()

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
    [f.close() for f in flist]
    # reset the env. variable
    os.environ['TBSLAS_RESULT_DIR'] = tbslas_result_dir

def run(
    ID,
    prog,
    pn_list,
    tl_list,
    dp_list,
    cq_list,
    ci_list,
    uf_list,
    np_list,
    nt_list,
    dt_list,
    tn_list,
    vs_list,
    mg_list,
    tt_list,
    di_list=None,
    ea_list=None):

    cmd_args = OrderedDict()
    cmd_args = generate_commands(prog,
                                 pn_list,
                                 tl_list,
                                 dp_list,
                                 cq_list,
                                 ci_list,
                                 uf_list,
                                 np_list,
                                 nt_list,
                                 dt_list,
                                 tn_list,
                                 vs_list,
                                 mg_list,
                                 tt_list,
                                 di_list,
                                 ea_list)
    execute_commands(cmd_args, ID)
