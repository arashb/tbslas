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
import pp

def generate_command_args(prog,\
                          pn_list,\
                          tl_list,\
                          dp_list,\
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
        ARGS    = ['-N'     , str(pn_list[counter]), \
                   '-tol'   , str(tl_list[counter]), \
                   '-d'     , str(dp_list[counter]), \
                   '-q'     , str(cq_list[counter]), \
                   '-cuf'   , str(uf_list[counter]), \
                   '-tn'    , str(2),              \
                   '-dt'    , str(1e-9),           \
                   '-vsr'   , str(0),                \
                   '-merge' , str(3),                \
                   '-omp'   , str(nt_list[counter])]
        if use_cubic:
            ARGS = ARGS + ['-cubic', '1']
        cmd_args[cmd_id] = utils.determine_command_prefix(np_list[counter]) + ['valgrind', '--tool=massif','--massif-out-file='+'valgrind-cmd-'+str(cmd_id)+'-%p.mem'] + [EXEC] + ARGS
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
    prog  = 'advection'
    tl_list = [\
            1e-2,\
            1e-4,\
            1e-7,\
            ]
    dp_list = [\
            # 6,\
            # 8,\
            # 10,\
            15,\
                ]
    cq_list = [\
            4,\
            6,\
            14,\
            ]
    np_list = [\
            1,\
            # 2,\
            # 4,\
            # 8,\
            # 16,\
            # 32,\
            ]
    use_cubic     = True
    max_np        = max(np_list)
    num_pnts      = 8**(math.floor(math.log(max_np,8)+1))
    uf            = 2
    table_counter = 0
    EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, prog)
    cmd_id = 0
    cmd_args = OrderedDict()
    for cq in cq_list:
        for tl in tl_list:
            if cq is 4 and tl is 1e-7:
                continue
            for dp in dp_list:
                # USE UF 4 FOR Q 14
                if cq is 14:
                    uf = 4
                for np in np_list:
                    ARGS    = ['-N'     , str(num_pnts), \
                               '-tol'   , str(tl), \
                               '-d'     , str(dp), \
                               '-q'     , str(cq), \
                               '-cuf'   , str(uf), \
                               '-tn'    , str(2),              \
                               '-dt'    , str(1e-9),           \
                               '-vs'    , str(1),                \
                               '-merge' , str(3),                \
                               '-omp'   , str(omp_num_threads)]
                    if use_cubic:
                        ARGS = ARGS + ['-cubic', '1']
                        cmd_args[cmd_id] = utils.determine_command_prefix(np) + ['valgrind', '--tool=massif','--massif-out-file='+'valgrind-cmd-'+str(cmd_id)+'-%p.mem'] + [EXEC] + ARGS
                    cmd_id = cmd_id + 1
    utils.execute_commands(cmd_args, prog+'-table-'+str(table_counter))
    table_counter = table_counter + 1


    raw_files_list = pp.list_raw_files(os.getcwd(), '.mem')
    print raw_files_list
    for raw_file in raw_files_list:
        os.popen('ms_print '+raw_file+' > '+raw_file+'.ms')
