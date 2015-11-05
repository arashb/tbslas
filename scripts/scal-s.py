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

TREE_MAX_DEPTH = 8

def generate_command_args(tl_list, \
                          dt_list, \
                          tn_list, \
                          np_list, \
                          merge_type,\
                          num_steps):

    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "advdiff-ss")

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS    = ['-N'    , str(8**5 ),              \
                   '-tol'  , str(tl_list[counter]),   \
                   '-d'    , str(TREE_MAX_DEPTH),     \
                   '-dt'   , str(2*dt_list[counter]), \
                   '-tn'   , str(tn_list[counter]),   \
                   '-vs'   , str(1),                  \
                   '-omp'  , str(nt_list[counter]),   \
                   '-cubic', str(1),                  \
                   '-cuf'  , str(4),                  \
                   '-merge', str(merge_type),         \
                   ]
        cmd_args[cmd_id] = utils.determine_command_prefix(np_list[counter]) + [EXEC] + ARGS
        cmd_id = cmd_id + 1
    return cmd_args

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    mpi_num_procs, omp_num_threads = utils.parse_args()

    num_steps = 1
    ############################################################################
    # TEST 1: STRONG SCALING
    ############################################################################
    tl_factor = 1
    tl_init   = 1e-0
    tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0, num_steps)]

    dt_factor = 1
    dt_init   = 0.25
    dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0, num_steps)]

    tn_factor = 1
    tn_init   = 1
    tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0, num_steps)]

    merge_type = 1

    # NUM MPI PROCESSES
    np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(tl_list, \
                                     dt_list, \
                                     tn_list, \
                                     np_list, \
                                     merge_type,\
                                     num_steps)
    utils.execute_commands(cmd_args,'merge-type-1')

    ############################################################################
    # TEST 2: STRONG SCALING
    ############################################################################
    tl_factor = 1
    tl_init   = 1e-0
    tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0, num_steps)]

    dt_factor = 1
    dt_init   = 0.25
    dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0, num_steps)]

    tn_factor = 1
    tn_init   = 1
    tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0, num_steps)]

    merge_type = 2

    # NUM MPI PROCESSES
    np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(tl_list, \
                                     dt_list, \
                                     tn_list, \
                                     np_list, \
                                     merge_type,\
                                     num_steps)

    utils.execute_commands(cmd_args,'merge-type-2')

    ############################################################################
    # TEST 3: STRONG SCALING
    ############################################################################
    tl_factor = 1
    tl_init   = 1e-0
    tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0, num_steps)]

    dt_factor = 1
    dt_init   = 0.25
    dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0, num_steps)]

    tn_factor = 1
    tn_init   = 1
    tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0, num_steps)]

    merge_type = 3

    # NUM MPI PROCESSES
    np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(tl_list, \
                                     dt_list, \
                                     tn_list, \
                                     np_list, \
                                     merge_type,\
                                     num_steps)

    utils.execute_commands(cmd_args,'merge-type-3')
