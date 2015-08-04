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

def generate_command_args(tl_init, tl_factor, \
                          dt_init, dt_factor, \
                          tn_init, tn_factor, \
                          np_list, np_factor, \
                          num_steps):
    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "advdiff-ss")

    tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0, num_steps)]
    dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0, num_steps)]
    tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0, num_steps)]
    # np_list = [np_init*math.pow(np_factor,float(cnt)) for cnt in range(0, num_steps)]
    np_list = [utils.MPI_TOTAL_NUM_PORCESSES          for cnt in range(0, num_steps)]
    nt_list = [utils.OMP_NUM_THREADS                  for cnt in range(0, num_steps)]

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS    = ['-N'   , str(8**5 ), \
                   '-tol' , str(tl_list[counter]),                              \
                   '-d', str(8), \
                   '-dt'  , str(2*dt_list[counter]),                            \
                   '-tn'  , str(tn_list[counter]),                              \
                   '-vs'  , str(1),                                             \
                   '-omp' , str(nt_list[counter]),
                   '-cubic', str(1),
                   '-cuf', str(4)]
        cmd_args[cmd_id] = utils.determine_command_prefix(np_list[counter]) + [EXEC] + ARGS
        cmd_id = cmd_id + 1
    return cmd_args

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    utils.prepare_environment(utils.OUTPUT_PREFIX)
    TOL_NUM_STEPS = 1
    if len(sys.argv) >= 4:
        TOL_NUM_STEPS   = int(sys.argv[3])

    ############################################################################
    # TEST 1: STRONG SCALING
    ############################################################################
    tl_factor = 1
    tl_init   = 1e-0
    dt_factor = 1
    dt_init   = 0.25
    tn_factor = 1
    tn_init   = 1
    np_factor = 2
    np_init   = 1

    cmd_args = generate_command_args(tl_init, tl_factor, \
                                     dt_init, dt_factor, \
                                     tn_init, tn_factor, \
                                     np_init, np_factor, \
                                     TOL_NUM_STEPS)
    utils.execute_commands(cmd_args,'table1')
