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

def generate_command_args(tl_list, \
                          dt_list, \
                          tn_list, \
                          de_list, \
                          q_list,  \
                          np_list, \
                          nt_list, \
                          num_steps):

    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "advection")

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS    = ['-N'   , str(8**math.ceil(math.log(np_list[counter],8))),    \
                   '-tol' , str(tl_list[counter]),                              \
                   '-q'   , str(q_list[counter]),                               \
                   '-dt'  , str(dt_list[counter]),                              \
                   '-tn'  , str(tn_list[counter]),                              \
                   '-d'   , str(de_list[counter]),                              \
                   '-omp' , str(nt_list[counter]),                              \
                   # '-vs'  , str(1),                               \
                   '-test', '2']
        cmd_args[cmd_id] = utils.determine_command_prefix(np_list[counter]) + [EXEC] + ARGS
        cmd_id = cmd_id + 1
    return cmd_args

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    mpi_num_procs, omp_num_threads = utils.parse_args()
    ############################################################################
    # TEST 1: TEMPORAL/SPATIAL ERROR
    ############################################################################
    num_steps = 6

    # TREE TOLERANCE
    tl_factor = 1
    tl_init   = 1e-10
    tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0, num_steps)]

    # TREE DEPTH
    # de_factor = 0
    # de_init   = 5
    # de_list = [de_init+cnt*de_factor                  for cnt in range(0, num_steps)]
    de_list = [5,5,5,7,7,7]

    # TIME RESOLUTION
    dt_factor = 1
    dt_init   = 6.28*1e-2
    dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0, num_steps)]

    # NUM TIME STEPS
    NUM_ROT = 5
    T_END = 6.28*NUM_ROT
    tn_factor = 1
    tn_init   = T_END/dt_init
    tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0, num_steps)]

    # CHEBYSHEV DEGREE
    # q_factor  = 1;
    # q_init    = 2;
    # q_list  = [q_init*math.pow(q_factor,float(cnt))   for cnt in range(0, num_steps)]
    q_list = [3,8,14,3,8,14]

    # NUM MPI PROCESSES
    np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(tl_list,\
                                     dt_list,\
                                     tn_list,\
                                     de_list,\
                                     q_list, \
                                     np_list,\
                                     nt_list,\
                                     num_steps)

    utils.execute_commands(cmd_args,'table1-ROT'+str(NUM_ROT))
