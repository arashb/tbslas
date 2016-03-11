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
import json
import utils

def scalws():
    ############################################################################
    # TEST 1: WEAK/STRONG SCALING FOR ADVECTION SOLVER
    ############################################################################
    prog  = 'advection'
    tl_list = [
            # 1e+0,
            1e-2,
            1e-4,
            1e-7,
            ]
    dp_list = [
            # 6,
            # 8,
            # 10,
            15,
        ]
    cq_list = [
            4,
            6,
            14,
            ]
    np_list = [
            1,
            2,
            4,
            8,
            16,
            32,
            ]
    # dt = 7.85398e-03
    # tn = 100
    dt = 0.0628
    tn = 200
    use_cubic     = True
    vsr = 0
    merge_type    = 3

    max_np        = max(np_list)
    num_pnts      = 8**(math.floor(math.log(max_np,8)+1))
    uf            = 2

    table_counter = 0
    for cq in cq_list:
        for tl in tl_list:
            if cq is 4 and tl is 1e-7:
                continue
            for dp in dp_list:
                # USE UF 4 FOR Q 14
                if cq is 14:
                    uf = 4
                cmd_args = OrderedDict()
                cmd_id = 0
                for np in np_list:
                    cmd_args[cmd_id] = utils.generate_commands(
                        prog    = prog,
                        pn_list = [num_pnts        ],
                        tl_list = [tl              ],
                        dp_list = [dp              ],
                        cq_list = [cq              ],
                        ci_list = [use_cubic       ],
                        uf_list = [uf              ],
                        np_list = [np              ],
                        nt_list = [omp_num_threads ],
                        dt_list = [dt              ],
                        tn_list = [tn              ],
                        vs_list = [vsr   ],
                        mg_list = [merge_type      ]
                        )[1]
                    cmd_id = cmd_id + 1
                utils.execute_commands(cmd_args, prog+'-table-'+str(table_counter))
                table_counter = table_counter + 1

def omp():
    ############################################################################
    # TEST: OMP SCALING
    ############################################################################
    mpi_num_procs, omp_num_threads = utils.parse_args()
    # prog     = 'advection'
    prog     = 'advdiff-ss'
    dt       = 0.0628
    vsr      = 0
    mrg_type = 3
    num_pnts = 8**(math.floor(math.log(mpi_num_procs,8)+1))

    num_steps = omp_num_threads
    nt_list = [cnt+1    for cnt in range(0,num_steps)]
    np_list = [mpi_num_procs for cnt in range(0,num_steps)]

    tn_list = [10       for cnt in range(0,num_steps)]
    pn_list = [num_pnts for cnt in range(0,num_steps)]
    tl_list = [1e-30    for cnt in range(0,num_steps)]
    dp_list = [4        for cnt in range(0,num_steps)]
    mg_list = [mrg_type for cnt in range(0,num_steps)]
    vs_list = [vsr      for cnt in range(0,num_steps)]
    cq_list = [14       for cnt in range(0,num_steps)]
    ci_list = [True     for cnt in range(0,num_steps)]
    uf_list = [4        for cnt in range(0,num_steps)]
    dt_list = [dt       for cnt in range(0,num_steps)]

    cmd_args = OrderedDict()
    cmd_args = utils.generate_commands(
        prog,
        pn_list,
        tl_list,
        dp_list,
        cq_list,
        ci_list,
        uf_list,
        np_list,
        nt_list,
        dt_list,
        tn_list,
        vs_list,
        mg_list)
    utils.execute_commands(cmd_args, 'omp')

def scals():
    mpi_num_procs, omp_num_threads = utils.parse_args()
    ############################################################################
    # TEST 1: STRONG SCALING
    ############################################################################
    # prog  = 'advection'
    prog  = 'advdiff-ss'
    np_list = [
             1,
             2,
             4,
             8,
             16,
             32,
             64,
            ]
    mrg_type = 3;
    max_np        = max(np_list)
    num_pnts      = 8**5#(math.floor(math.log(max_np,8)+1))

    num_steps = len(np_list)
    pn_list = [num_pnts        for cnt in range(0,num_steps)]
    tl_list = [1e-0            for cnt in range(0,num_steps)]
    dp_list = [15              for cnt in range(0,num_steps)]
    cq_list = [14              for cnt in range(0,num_steps)]
    ci_list = [True            for cnt in range(0,num_steps)]
    uf_list = [4               for cnt in range(0,num_steps)]
    nt_list = [omp_num_threads for cnt in range(0,num_steps)]
    dt_list = [0.125           for cnt in range(0,num_steps)]
    tn_list = [10              for cnt in range(0,num_steps)]
    vs_list = [0               for cnt in range(0,num_steps)]
    mg_list = [mrg_type        for cnt in range(0,num_steps)]

    cmd_args = OrderedDict()
    cmd_args = utils.generate_commands(
        prog,
        pn_list,
        tl_list,
        dp_list,
        cq_list,
        ci_list,
        uf_list,
        np_list,
        nt_list,
        dt_list,
        tn_list,
        vs_list,
        mg_list)
    utils.execute_commands(cmd_args, 'sscal')

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    # omp()
    # scalws()
    scals()
