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
                          de_init, de_factor, \
                          q_init , q_factor,  \
                          num_steps):
    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "advection")

    tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0, num_steps)]
    dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0, num_steps)]
    tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0, num_steps)]
    de_list = [de_init+cnt*de_factor                  for cnt in range(0, num_steps)]
    q_list  = [q_init*math.pow(q_factor,float(cnt))   for cnt in range(0, num_steps)]
    np_list = [utils.MPI_TOTAL_NUM_PORCESSES          for cnt in range(0, num_steps)]
    nt_list = [utils.OMP_NUM_THREADS                  for cnt in range(0, num_steps)]

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
    utils.prepare_environment(utils.OUTPUT_PREFIX)
    ############################################################################
    # TEST 1: SPATIAL ERROR
    ############################################################################
    # TOL_NUM_STEPS = 7
    # if len(sys.argv) >= 4:
    #     TOL_NUM_STEPS   = int(sys.argv[3])
    # tl_factor = 0.1
    # tl_init   = 1e-10
    # de_factor = 1
    # de_init   = 4;
    # dt_factor = 1
    # dt_init   = 1e-3
    # tn_factor = 1.0/dt_factor
    # tn_init   = 1
    # q_factor  = 1;
    # q_init    = 14;
    # cmd_args = generate_command_args(tl_init, tl_factor, \
    #                                  dt_init, dt_factor, \
    #                                  tn_init, tn_factor, \
    #                                  de_init, de_factor, \
    #                                  q_init, q_factor, \
    #                                  TOL_NUM_STEPS)
    # utils.execute_commands(cmd_args,'table1')
    ############################################################################
    # TEST 2: TEMPORAL/SPATIAL ERROR
    ############################################################################
    TOL_NUM_STEPS = 1
    NUM_ROT = 10
    T_END = 6.28*NUM_ROT
    tl_factor = 1
    tl_init   = 1e-10
    de_factor = 0
    de_init   = 6
    dt_factor = 1
    dt_init   = 6.28*1e-2
    tn_factor = 1
    tn_init   = T_END/dt_init
    q_factor  = 1;
    q_init    = 2;
    cmd_args = generate_command_args(tl_init, tl_factor, \
                                         dt_init, dt_factor, \
                                         tn_init, tn_factor, \
                                         de_init, de_factor, \
                                         q_init, q_factor, \
                                         TOL_NUM_STEPS)
    utils.execute_commands(cmd_args,'table2-ROT'+str(NUM_ROT))
