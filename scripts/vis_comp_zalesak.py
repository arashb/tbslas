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

def generate_command_args(max_depth, mpi_num_procs, omp_num_threads, use_cubic):

    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "advection")

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1
    ARGS    = ['-N'   , '4096', \
               '-tol' , '1e-2', \
               '-d'   , str(max_depth), \
               '-dt'  , '0.0628', \
               '-tn'  , '500', \
               '-test', str(2), \
               '-omp' , str(omp_num_threads), \
               '-q'   , str(4), \
               # '-vs'  , '1'
               ]
    if use_cubic:
        ARGS = ARGS + ['-cubic', '1']
        ARGS = ARGS + ['-cuf', '2']

    cmd_args[cmd_id] = utils.determine_command_prefix(mpi_num_procs) + [EXEC] + ARGS

    return cmd_args

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    mpi_num_procs, omp_num_threads = utils.parse_args()
    max_depth = 8
    ############################################################################
    # TEST 1: TEMPORAL CONVERGENCE
    ############################################################################
    use_cubic = True
    cmd_args = generate_command_args(max_depth, mpi_num_procs, omp_num_threads, use_cubic)
    utils.execute_commands(cmd_args, 'vis-zalesak')
