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
import json

def generate_command_args(prog,
                          tl_list,\
                          cq_list,\
                          uf_list,\
                          use_cubic,\
                          np_list,\
                          nt_list,\
                          num_steps):

    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, prog)

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS    = ['-N'     , str(8), \
                   '-tol'   , str(tl_list[counter]),                              \
                   '-q'     , str(cq_list[counter]),                              \
                   '-cuf'   , str(uf_list[counter]),                              \
                   '-tn'    , str(500),                              \
                   '-dt'    , str(0.0628),                              \
                   '-vs'  , str(1),                               \
                   '-merge' , str(3),\
                   '-omp'   , str(nt_list[counter])]
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
    tl_list = [\
            1e-2,\
            1e-4,\
            1e-7,\
            ]
    cq_list = [\
            4,\
            6,\
            14,\
            ]
    np_list = [\
            1,\
            16,\
            32,\
            ]

    use_cubic = True

    cmd_args = OrderedDict()
    cmd_id = 0
    uf = 2
    for cq in cq_list:
        for tl in tl_list:
            if cq is 14:
                uf = 4
            # if cq is 4 and tl is 1e-7:
            #     continue
            for np in np_list:
                cmd_args[cmd_id] = generate_command_args(\
                                                            prog      = 'advection',\
                                                            tl_list   = [tl],\
                                                            cq_list   = [cq],\
                                                            uf_list   = [uf],\
                                                            use_cubic = True,\
                                                            np_list   = [np],\
                                                            nt_list   = [omp_num_threads],\
                                                            num_steps = 1)[1]
                cmd_id = cmd_id + 1

    # print(json.dumps(cmd_args, indent=4))
    utils.execute_commands(cmd_args, 'table1-advection')

    ############################################################################
    # TEST 2:
    ############################################################################
    tl_list = [\
            1e-2,\
            1e-4,\
            1e-7,\
            ]
    cq_list = [\
            4,\
            6,\
            14,\
            ]

    np_list = [\
               1,\
               16,\
               32,\
                ]

    use_cubic = True

    cmd_args = OrderedDict()
    cmd_id = 0
    uf = 2
    for cq in cq_list:
        for tl in tl_list:
            if cq is 14:
                uf = 4
            # if cq is 4 and tl is 1e-7:
            #     continue
            for np in np_list:
                cmd_args[cmd_id] = generate_command_args(\
                                                            prog      = 'advdiff-ss',\
                                                            tl_list   = [tl],\
                                                            cq_list   = [cq],\
                                                            uf_list   = [uf],\
                                                            use_cubic = True,\
                                                            np_list   = [np],\
                                                            nt_list   = [omp_num_threads],\
                                                            num_steps = 1)[1]
                cmd_id = cmd_id + 1

    # print(json.dumps(cmd_args, indent=4))
    utils.execute_commands(cmd_args, 'table2-advdiff')
