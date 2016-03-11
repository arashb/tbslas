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
                          merge_type, \
                          np_list, \
                          nt_list, \
                          num_steps):

    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "advection")

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS    = ['-N'   , str(8**3), \
                   '-tol' , str(1e-10),                                      \
                   '-q'   , str(5),                                          \
                   '-d'   , str(de_list[counter]),                           \
                   '-dt'  , str(dt_list[counter]),                           \
                   '-tn'  , str(tn_list[counter]),                           \
                   '-test', str(test_list[counter]),                         \
                   '-merge', str(merge_type),                                \
                   '-vsr'  , str(0),                                          \
                   # '-cubic',str(1),                                        \
                   # '-cuf'  ,str(8),                                        \
                   '-omp' , str(nt_list[counter])]
        cmd_args[cmd_id] = utils.determine_command_prefix(np_list[counter]) + [EXEC] + ARGS
        cmd_id = cmd_id + 1
    return cmd_args

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    mpi_num_procs, omp_num_threads = utils.parse_args()
    for merge_type in range(1,4):
        # ######################################################################
        # # TEST 3: V,C depth: [5, 7, 9] config: irregular V, irregular C
        # ######################################################################
        np_list = [
                1,
                2,
                4,
                8,
                16,
                32,
                ]

        num_steps = len(np_list)
        de_factor = 1
        de_init   = 6
        de_list   = [de_init for cnt in range(0,num_steps)]

        nt_list   = [omp_num_threads  for cnt in range(0,num_steps)]

        dt_init   = 0.25
        dt_list   = [dt_init for cnt in range(0,num_steps)]

        tn_init   = 1
        tn_list   = [tn_init for cnt in range(0,num_steps)]

        test_init = 6
        test_list = [test_init for cnt in range(0,num_steps)]

        cmd_args = generate_command_args(de_list,
                                         dt_list,
                                         tn_list,
                                         test_list,
                                         merge_type,
                                         np_list,
                                         nt_list,
                                         num_steps)
        utils.execute_commands(cmd_args, 'test-3-merge-type-'+str(merge_type))
