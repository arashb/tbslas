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
                          cq_list,\
                          uf_list,\
                          use_cubic,\
                          np_list,\
                          nt_list,\
                          num_steps):

    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "advection")

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS    = ['-N'   , str(8), \
                   '-tol' , str(tl_list[counter]),                              \
                   '-q'   , str(cq_list[counter]),                              \
                   '-cuf' , str(uf_list[counter]),                              \
                   '-tn'  , str(10),                              \
                   '-dt'  , str(0.015625),                              \
                   # '-vs'  , str(1),                               \
                   '-omp' , str(nt_list[counter])]
        if use_cubic:
            ARGS = ARGS + ['-cubic', '1']
        cmd_args[cmd_id] = utils.determine_command_prefix(np_list[counter]) + [EXEC] + ARGS
        cmd_id = cmd_id + 1
    return cmd_args

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    mpi_num_procs, omp_num_threads = utils.parse_args()

    ############################################################################
    # TEST 1:
    ############################################################################
    tl_list = [1e-2, 1e-4, 1e-7]
    cq_list = [6 , 10, 14]
    uf_list = [2, 4]

    use_cubic = True

    cmd_args = OrderedDict()
    cmd_id = 0
    for tl in tl_list:
        for cq in cq_list:
            if cq is 3:
                my_uf_list = [1]
            else:
                my_uf_list = uf_list
            for uf in my_uf_list:
                cmd_args[cmd_id] = generate_command_args([tl],\
                                                 [cq],\
                                                 [uf],\
                                                 use_cubic,\
                                                 [mpi_num_procs],\
                                                 [omp_num_threads],\
                                                1)[1]
                cmd_id = cmd_id + 1
    print cmd_args
    utils.execute_commands(cmd_args, 'table1')

    # ##########################################################################
    # # TEST 2:
    # ##########################################################################
    tl_list = [1e-2, 1e-4, 1e-7]
    cq_list = [3, 6 , 10, 14]
    uf_list = [1]

    use_cubic = False

    cmd_args = OrderedDict()
    cmd_id = 0
    for tl in tl_list:
        for cq in cq_list:
            if cq is 3:
                my_uf_list = [1]
            else:
                my_uf_list = uf_list
            for uf in my_uf_list:
                cmd_args[cmd_id] = generate_command_args([tl],\
                                                 [cq],\
                                                 [uf],\
                                                 use_cubic,\
                                                 [mpi_num_procs],\
                                                 [omp_num_threads],\
                                                1)[1]
                cmd_id = cmd_id + 1
    print cmd_args
    utils.execute_commands(cmd_args, 'table2')
