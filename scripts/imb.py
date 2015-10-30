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
        ARGS    = ['-N'   , str(8**math.ceil(math.log(np_list[counter],8))), \
                   '-tol' , str(1e-10),                                      \
                   '-q'   , str(5),                                          \
                   '-d'   , str(de_list[counter]),                           \
                   '-dt'  , str(dt_list[counter]),                           \
                   '-tn'  , str(tn_list[counter]),                           \
                   '-test', str(test_list[counter]),                         \
                   '-merge', str(merge_type),                                \
                   '-vs'  , str(1),                                          \
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
    num_steps = 2
    for merge_type in range(1,4):
    # ############################################################################
    # # TEST 1: V,C depth: [6] config: regular V, regular C
    # ############################################################################
        # de_factor = 1
        # de_init   = 5
        # dt_factor = 1
        # dt_init   = 0.25
        # tn_factor = 1.0
        # tn_init   = 1
        # test_init = 4
        # test_factor = 0;
        # cmd_args = generate_command_args(de_init, de_factor, \
        #                              dt_init, dt_factor, \
        #                              tn_init, tn_factor, \
        #                              test_init, test_factor, \
        #                              merge_type, \
        #                              1)
        # utils.execute_commands(cmd_args, 'test-1-merge-type-'+str(merge_type))

    # ############################################################################
    # # TEST 2: V depth: [6] C depth: [5, 7, 9] config: regular V, irregular C
    # ############################################################################
        # de_factor = 1
        # de_init   = 5
        # dt_factor = 1
        # dt_init   = 0.25
        # tn_factor = 1.0
        # tn_init   = 1
        # test_init = 5
        # test_factor = 0;
        # cmd_args = generate_command_args(de_init, de_factor, \
        #                              dt_init, dt_factor, \
        #                              tn_init, tn_factor, \
        #                              test_init, test_factor, \
        #                              merge_type, \
        #                              num_steps)
        # utils.execute_commands(cmd_args, 'test-2-merge-type-'+str(merge_type))

    ############################################################################
    # TEST 3: V,C depth: [5, 7, 9] config: irregular V, irregular C
    ############################################################################
        de_factor = 1
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

        cmd_args = generate_command_args(de_list, \
                                        dt_list, \
                                        tn_list, \
                                        test_list, \
                                        merge_type, \
                                        np_list, \
                                        nt_list, \
                                        num_steps)

        utils.execute_commands(cmd_args, 'test-3-merge-type-'+str(merge_type))
