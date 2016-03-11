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
                          # de_list, \
                          # q_list,  \
                          np_list, \
                          nt_list, \
                          num_steps):

    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "traj")

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS    = ['-N'   , str(8**math.ceil(math.log(np_list[counter],8))), \
                   '-tol' , str(tl_list[counter]),                           \
                   '-dt'  , str(dt_list[counter]),                           \
                   '-tn'  , str(tn_list[counter]),                           \
                   '-vsr' , str(0),                                          \
                   '-omp' , str(nt_list[counter])]
        cmd_args[cmd_id] = utils.determine_command_prefix(np_list[counter]) + [EXEC] + ARGS
        cmd_id = cmd_id + 1
    return cmd_args

def test1():
    mpi_num_procs, omp_num_threads = utils.parse_args()
    num_steps = 8
    T_END     = 2*math.pi
    ############################################################################
    # TEST 1: TEMPORAL ERROR
    ############################################################################
    tl_fact = 1#0.1
    tl_init = 1e-5
    tl_list = [tl_init*math.pow(tl_fact,float(cnt)) for cnt in range(0,num_steps)]

    tn_fact = 2;
    tn_init = 50
    tn_list = [tn_init*math.pow(tn_fact,float(cnt)) for cnt in range(0,num_steps)]

    dt_fact = 0.5
    dt_init = T_END/tn_init
    dt_list = [dt_init*math.pow(dt_fact,float(cnt)) for cnt in range(0,num_steps)]


    # NUM MPI PROCESSES
    np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(tl_list,\
                                     dt_list,\
                                     tn_list,\
                                     # de_list,\
                                     # q_list, \
                                     np_list,\
                                     nt_list,\
                                     num_steps)

    utils.execute_commands(cmd_args, 'table1')

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    test1()
