#!/bin/env python
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

def generate_command_args(tl_init, tl_factor, cuf_init, cuf_factor, use_cubic, \
                          use_anal, num_steps):
    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "cubic")

    tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0,num_steps)]
    uf_list = [cuf_init*math.pow(cuf_factor,float(cnt)) for cnt in range(0,num_steps)]
    np_list = [utils.MPI_TOTAL_NUM_PORCESSES          for cnt in range(0,TOL_NUM_STEPS)]
    nt_list = [utils.OMP_NUM_THREADS                  for cnt in range(0,TOL_NUM_STEPS)]

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS    = ['-N'   , str(8**math.ceil(math.log(np_list[counter],8))), \
                   '-tol' , str(tl_list[counter]),                              \
                   '-cuf' , str(uf_list[counter]),                              \
                   '-vs'  , str(1),                               \
                   '-omp' , str(nt_list[counter])]
        if use_cubic:
            ARGS = ARGS + ['-cubic', '1']
        if use_anal:
            ARGS = ARGS + ['-ca','1']
        cmd_args[cmd_id] = utils.determine_command_prefix(np_list[counter]) + [EXEC] + ARGS
        cmd_id = cmd_id + 1
    return cmd_args

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    utils.prepare_environment(utils.OUTPUT_PREFIX)
    TOL_NUM_STEPS = 8
    if len(sys.argv) >= 4:
        TOL_NUM_STEPS   = int(sys.argv[3])
    ############################################################################
    # TEST 1:
    ############################################################################
    tl_factor  = 1
    tl_init    = 1e-8
    cuf_factor = 2
    cuf_init   = 2
    use_cubic  = True
    use_anal   = True
    cmd_args = generate_command_args(tl_init, tl_factor, \
                                     cuf_init, cuf_factor, \
                                     use_cubic, use_anal, \
                                     TOL_NUM_STEPS)
    utils.execute_commands(cmd_args, 'table1')
    ############################################################################
    # TEST 2:
    ############################################################################
    tl_factor  = 0.1
    tl_init    = 1e-6
    cuf_factor = 1
    cuf_init   = 2
    use_cubic  = False
    use_anal   = False
    cmd_args = generate_command_args(tl_init, tl_factor, \
                                     cuf_init, cuf_factor, \
                                     use_cubic, use_anal, \
                                     TOL_NUM_STEPS)
    utils.execute_commands(cmd_args, 'table2')
    ############################################################################
    # TEST 3:
    ############################################################################
    tl_factor  = 0.1
    tl_init    = 1e-3
    cuf_factor = 2
    cuf_init   = 2
    use_cubic  = True
    use_anal   = False
    cmd_args = generate_command_args(tl_init, tl_factor, \
                                     cuf_init, cuf_factor, \
                                     use_cubic, use_anal, \
                                     TOL_NUM_STEPS)
    utils.execute_commands(cmd_args, 'table3')
