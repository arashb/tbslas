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


def generate_command_args(de_list, \
                          dt_list, \
                          tn_list, \
                          test_list, \
                          mt_list, \
                          np_list, \
                          nt_list, \
                          num_steps):

    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "advection")

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS    = ['-N'   , str(8**math.ceil(math.log(np_list[counter],8))), \
                   '-tol' , str(1e-10),                                      \
                   '-q'   , str(5),                                          \
                   '-d'   , str(de_list[counter]),                           \
                   '-dt'  , str(dt_list[counter]),                           \
                   '-tn'  , str(tn_list[counter]),                           \
                   '-test', str(test_list[counter]),                         \
                   '-merge', str(mt_list[counter]),                                \
                   '-vsr'  , str(0),                                          \
                   # '-cubic',str(1),                                        \
                   # '-cuf'  ,str(8),                                        \
                   '-omp' , str(nt_list[counter])]
        cmd_args[cmd_id] = utils.determine_command_prefix(np_list[counter]) + [EXEC] + ARGS
        cmd_id = cmd_id + 1
    return cmd_args

def test1():
    ############################################################################
    # TEST 1: V,C depth: [6] config: regular V, regular C
    ############################################################################
    mpi_num_procs, omp_num_threads = utils.parse_args()
    num_steps = 3

    mt_factor = 1
    mt_init = 1
    mt_list = [mt_init+cnt*mt_factor                  for cnt in range(0, num_steps)]

    de_factor = 0
    de_init   = 5
    de_list = [de_init+cnt*de_factor                  for cnt in range(0, num_steps)]

    dt_factor = 1
    dt_init   = 0.25
    dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0,num_steps)]

    tn_factor = 1.0
    tn_init   = 1
    tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0,num_steps)]

    test_init = 4
    test_factor = 0;
    test_list = [test_init+cnt*test_factor            for cnt in range(0,num_steps)]

    # NUM MPI PROCESSES
    np_list = [mpi_num_procs  for cnt in range(0, num_steps)]
    # NUM OMP THREADS
    nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(de_list,
                                     dt_list,
                                     tn_list,
                                     test_list,
                                     mt_list,
                                     np_list,
                                     nt_list,
                                     num_steps)
    utils.execute_commands(cmd_args, 'test1')

def test2():
    ############################################################################
    # TEST 2: V depth: [6] C depth: [5, 7, 9] config: regular V, irregular C
    ############################################################################
    mpi_num_procs, omp_num_threads = utils.parse_args()
    num_steps = 3

    mt_factor = 1
    mt_init = 1
    mt_list = [mt_init+cnt*mt_factor                  for cnt in range(0, num_steps)]

    de_factor = 0
    de_init   = 5
    de_list = [de_init+cnt*de_factor                  for cnt in range(0, num_steps)]

    dt_factor = 1
    dt_init   = 0.25
    dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0,num_steps)]

    tn_factor = 1.0
    tn_init   = 1
    tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0,num_steps)]

    test_init = 5
    test_factor = 0;
    test_list = [test_init+cnt*test_factor            for cnt in range(0,num_steps)]

    # NUM MPI PROCESSES
    np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(de_list,
                                     dt_list,
                                     tn_list,
                                     test_list,
                                     mt_list,
                                     np_list,
                                     nt_list,
                                     num_steps)
    utils.execute_commands(cmd_args, 'test2')

def test3():
    ############################################################################
    # TEST 3: V,C depth: [5, 7, 9] config: irregular V, irregular C
    ############################################################################
    mpi_num_procs, omp_num_threads = utils.parse_args()
    num_steps = 3

    mt_factor = 1
    mt_init = 1
    mt_list = [mt_init+cnt*mt_factor                  for cnt in range(0, num_steps)]

    de_factor = 0
    de_init   = 5
    de_list = [de_init+cnt*de_factor                  for cnt in range(0, num_steps)]

    dt_factor = 1
    dt_init   = 0.25
    dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0,num_steps)]

    tn_factor = 1.0
    tn_init   = 1
    tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0,num_steps)]

    test_init = 6
    test_factor = 0;
    test_list = [test_init+cnt*test_factor            for cnt in range(0,num_steps)]

    # NUM MPI PROCESSES
    np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(de_list,
                                     dt_list,
                                     tn_list,
                                     test_list,
                                     mt_list,
                                     np_list,
                                     nt_list,
                                     num_steps)
    utils.execute_commands(cmd_args, 'test3')

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    # test1()
    # test2()
    test3()
