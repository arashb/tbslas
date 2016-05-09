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
import os
import subprocess
import math
import sys
from collections import OrderedDict
import utils

def generate_command_args(tl_list,\
                          uf_list,\
                          use_cubic,\
                          np_list,\
                          nt_list,\
                          num_steps):

    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "field-set")

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS    = ['-N'   , str(8**math.ceil(math.log(np_list[counter],8))), \
                   '-tol' , str(tl_list[counter]),                              \
                   '-cuf' , str(uf_list[counter]),                              \
                   '-vsr' , str(0),                               \
                   '-omp' , str(nt_list[counter])]
        if use_cubic:
            ARGS = ARGS + ['-cubic', '1']
        cmd_args[cmd_id] = utils.determine_command_prefix(np_list[counter]) + [EXEC] + ARGS
        cmd_id = cmd_id + 1
    return cmd_args

def test1():
    mpi_num_procs, omp_num_threads = utils.parse_args()
    num_steps = 8
    ############################################################################
    # TEST 1:
    ############################################################################
    tl_fact   = 1
    tl_init   = 1e-8
    tl_list   = [tl_init*math.pow(tl_fact,float(cnt)) for cnt in range(0,num_steps)]

    cuf_fact  = 2
    cuf_init  = 2
    uf_list   = [cuf_init*math.pow(cuf_fact,float(cnt)) for cnt in range(0,num_steps)]

    use_cubic = True

    # NUM MPI PROCESSES
    np_list   = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list   = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(tl_list,\
                                     uf_list,\
                                     use_cubic,\
                                     np_list,\
                                     nt_list,\
                                     num_steps)

    utils.execute_commands(cmd_args, 'table1')


def test2():
    mpi_num_procs, omp_num_threads = utils.parse_args()
    num_steps = 8
    ############################################################################
    # TEST 2:
    ############################################################################
    tl_fact   = 0.1
    tl_init   = 1e-6
    tl_list   = [tl_init*math.pow(tl_fact,float(cnt)) for cnt in range(0,num_steps)]

    cuf_fact  = 1
    cuf_init  = 2
    uf_list   = [cuf_init*math.pow(cuf_fact,float(cnt)) for cnt in range(0,num_steps)]

    use_cubic = False

    # NUM MPI PROCESSES
    np_list   = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list   = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(tl_list,\
                                     uf_list,\
                                     use_cubic,\
                                     np_list,\
                                     nt_list,\
                                     num_steps)

    utils.execute_commands(cmd_args, 'table2')

def test3():
    mpi_num_procs, omp_num_threads = utils.parse_args()
    num_steps = 8
    ############################################################################
    # TEST 3:
    ############################################################################
    tl_fact   = 0.1
    tl_init   = 1e-3
    tl_list   = [tl_init*math.pow(tl_fact,float(cnt)) for cnt in range(0,num_steps)]

    cuf_fact  = 2
    cuf_init  = 2
    uf_list   = [cuf_init*math.pow(cuf_fact,float(cnt)) for cnt in range(0,num_steps)]

    use_cubic = True
    # use_anal   = False

    # NUM MPI PROCESSES
    np_list   = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list   = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(tl_list,\
                                     uf_list,\
                                     use_cubic,\
                                     np_list,\
                                     nt_list,\
                                     num_steps)

    utils.execute_commands(cmd_args, 'table3')

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    tbslas_dir = os.environ['TBSLAS_RESULT_DIR']
    import time
    TIMESTR       = time.strftime("%Y%m%d-%H%M%S-")+str(time.time())
    os.environ['TBSLAS_RESULT_DIR'] =  os.path.join(tbslas_dir,'conv-cubic-'+TIMESTR)
    if not os.path.exists(os.environ['TBSLAS_RESULT_DIR']):
        os.makedirs(os.environ['TBSLAS_RESULT_DIR'])

    test1()
    test2()
    test3()

    os.environ['TBSLAS_RESULT_DIR'] = tbslas_dir
