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

TREE_MAX_DEPTH = 8

# def generate_command_args(tl_list, \
#                           dt_list, \
#                           tn_list, \
#                           np_list, \
#                           merge_type,\
#                           num_steps):

def generate_command_args(prog,\
                          pn_list,\
                          tl_list,\
                          cq_list,\
                          uf_list,\
                          dt_list,\
                          tn_list,\
                          use_cubic,\
                          merge_type,\
                          np_list,\
                          nt_list,\
                          num_steps):

    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, prog)

    # generate a dictionary data type of commands
    cmd_args = OrderedDict()
    cmd_id = 1;
    for counter in range(0,num_steps):
        ARGS    = ['-N'    , str(pn_list[counter] ), \
                   '-tol'  , str(tl_list[counter]),  \
                   '-q'    , str(cq_list[counter]),  \
                   '-d'    , str(TREE_MAX_DEPTH),    \
                   '-dt'   , str(dt_list[counter]),  \
                   '-tn'   , str(tn_list[counter]),  \
                   '-vs'   , str(1),                 \
                   '-omp'  , str(nt_list[counter]),  \
                   '-cuf'  , str(uf_list[counter]),  \
                   '-merge', str(merge_type),        \
                   ]
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
    # TEST 1: STRONG SCALING
    ############################################################################
    prog  = 'advection'
    tl_list = [1e-0]
    dt_list = [0.25]
    tn_list = [100]
    cq_list = [\
            # 4,\
            # 6,\
            14,\
            ]
    uf_list = [4]

    np_list = [\
            1,\
            2,\
            4,\
            8,\
            16,\
            # 32,\
            ]
    mt_list = [1,2,3]
    use_cubic     = True
    max_np        = max(np_list)
    num_pnts      = 8**(math.floor(math.log(max_np,8)+1))
    table_counter = 0
    for merge_type in mt_list:
            cmd_args = OrderedDict()
            cmd_id = 0
            for np in np_list:
                cmd_args[cmd_id] = generate_command_args(\
                                                         prog      = prog,\
                                                         pn_list   = [num_pnts],\
                                                         tl_list   = tl_list,\
                                                         cq_list   = cq_list,\
                                                         uf_list   = uf_list,\
                                                         dt_list   = dt_list,\
                                                         tn_list   = tn_list,\
                                                         use_cubic = True,\
                                                         merge_type= merge_type,\
                                                         np_list   = [np],\
                                                         nt_list   = [omp_num_threads],\
                                                         num_steps = 1)[1]
                cmd_id = cmd_id + 1
            # print(json.dumps(cmd_args, indent=4))
            utils.execute_commands(cmd_args, prog+'-table-'+str(table_counter))
            table_counter = table_counter + 1

# def generate_command_args(prog,\
#                           pn_list,\
#                           tl_list,\
#                           cq_list,\
#                           uf_list,\
#                           dt_list,\
#                           tn_list,\
#                           use_cubic,\
#                           merge_type,\
#                           np_list,\
#                           nt_list,\
#                           num_steps):

    # tl_factor = 1
    # tl_init   = 1e-0
    # tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0, num_steps)]

    # dt_factor = 1
    # dt_init   = 0.25
    # dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0, num_steps)]

    # tn_factor = 1
    # tn_init   = 1
    # tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0, num_steps)]

    # merge_type = 1

    # # NUM MPI PROCESSES
    # np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # # NUM OMP THREADS
    # nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    # cmd_args = generate_command_args(tl_list, \
    #                                  dt_list, \
    #                                  tn_list, \
    #                                  np_list, \
    #                                  merge_type,\
    #                                  num_steps)
    # utils.execute_commands(cmd_args,'merge-type-1')

    # ############################################################################
    # # TEST 2: STRONG SCALING
    # ############################################################################
    # tl_factor = 1
    # tl_init   = 1e-0
    # tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0, num_steps)]

    # dt_factor = 1
    # dt_init   = 0.25
    # dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0, num_steps)]

    # tn_factor = 1
    # tn_init   = 1
    # tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0, num_steps)]

    # merge_type = 2

    # # NUM MPI PROCESSES
    # np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # # NUM OMP THREADS
    # nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    # cmd_args = generate_command_args(tl_list, \
    #                                  dt_list, \
    #                                  tn_list, \
    #                                  np_list, \
    #                                  merge_type,\
    #                                  num_steps)

    # utils.execute_commands(cmd_args,'merge-type-2')

    # ############################################################################
    # # TEST 3: STRONG SCALING
    # ############################################################################
    # tl_factor = 1
    # tl_init   = 1e-0
    # tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0, num_steps)]

    # dt_factor = 1
    # dt_init   = 0.25
    # dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0, num_steps)]

    # tn_factor = 1
    # tn_init   = 1
    # tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0, num_steps)]

    # merge_type = 3

    # # NUM MPI PROCESSES
    # np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # # NUM OMP THREADS
    # nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    # cmd_args = generate_command_args(tl_list, \
    #                                  dt_list, \
    #                                  tn_list, \
    #                                  np_list, \
    #                                  merge_type,\
    #                                  num_steps)

    # utils.execute_commands(cmd_args,'merge-type-3')
