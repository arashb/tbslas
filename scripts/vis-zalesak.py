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
from utils import *

def generate_command_args(max_depth):
    EXEC = os.path.join(TBSLAS_EXAMPLES_BIN_DIR, "advection")
    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1
    ARGS    = ['-N'   , '4096', \
               '-tol' , '1e-3', \
               '-d'   , str(max_depth), \
               '-dt'  , '0.0628', \
               '-tn'  , '500', \
               '-test', str(2), \
               '-omp' , str(OMP_NUM_THREADS), \
               '-q'   , str(8), \
               # '-vs'  , '1'
               ]
    cmd_args[cmd_id] = [EXEC] + ARGS
    return cmd_args
################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    prepare_environment(OUTPUT_PREFIX)
    max_depth = 7
    if len(sys.argv) >= 4:
        max_depth   = int(sys.argv[3])
    ############################################################################
    # TEST 1: TEMPORAL CONVERGENCE
    ############################################################################
    cmd_args = generate_command_args(max_depth)
    execute_commands(cmd_args, 'vis-zalesak')
