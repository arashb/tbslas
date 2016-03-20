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

def sscal():
    ############################################################################
    # STRONG SCALING
    ############################################################################
    mpi_num_procs, omp_num_threads = utils.parse_args()
    for merge_type in range(1,4):
        prog  = 'advection'
        # prog  = 'advdiff-ss'
        np_list = [
             1,
             2,
             4,
             8,
             16,
             32,
             64,
            ]
        max_np        = max(np_list)
        num_pnts      = 8**3#(math.floor(math.log(max_np,8)+1))

        num_steps = len(np_list)
        pn_list = [num_pnts        for cnt in range(0,num_steps)]
        tl_list = [1e-5            for cnt in range(0,num_steps)]
        dp_list = [6               for cnt in range(0,num_steps)]
        cq_list = [5               for cnt in range(0,num_steps)]
        ci_list = [False           for cnt in range(0,num_steps)]
        uf_list = [4               for cnt in range(0,num_steps)]
        nt_list = [omp_num_threads for cnt in range(0,num_steps)]
        dt_list = [0.25            for cnt in range(0,num_steps)]
        tn_list = [1               for cnt in range(0,num_steps)]
        vs_list = [0               for cnt in range(0,num_steps)]
        mg_list = [merge_type      for cnt in range(0,num_steps)]
        tt_list = [6               for cnt in range(0,num_steps)]

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
            mg_list,
            tt_list)
        utils.execute_commands(cmd_args,  prog+'-sscal-mt-'+str(merge_type))

def wscal():
    mpi_num_procs, omp_num_threads = utils.parse_args()
    for merge_type in range(1,4):
        prog  = 'advection'
        # prog  = 'advdiff-ss'
        np_list = [
            1,
            2,
            4,
            8,
            16,
            32,
            64,
        ]
        max_np        = max(np_list)
        num_pnts      = 8**3#(math.floor(math.log(max_np,8)+1))

        num_steps = len(np_list)
        pn_list = [num_pnts        for cnt in range(0,num_steps)]

        dp_factor = 1
        dp_init   = 6
        dp_list   = [dp_init+cnt*dp_factor for cnt in range(0,num_steps)]

        tl_factor = 0.1
        tl_init   = 1e-5
        tl_list   = [tl_init*math.pow(tl_factor,float(cnt))     for cnt in range(0,num_steps)]

        cq_list = [5               for cnt in range(0,num_steps)]
        ci_list = [False           for cnt in range(0,num_steps)]
        uf_list = [4               for cnt in range(0,num_steps)]
        nt_list = [omp_num_threads for cnt in range(0,num_steps)]
        dt_list = [0.25            for cnt in range(0,num_steps)]
        tn_list = [1               for cnt in range(0,num_steps)]
        vs_list = [0               for cnt in range(0,num_steps)]
        mg_list = [merge_type      for cnt in range(0,num_steps)]
        tt_list = [6               for cnt in range(0,num_steps)]


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
            mg_list,
            tt_list)
        utils.execute_commands(cmd_args,  prog+'-wscal-mt-'+str(merge_type))

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    sscal()
    wscal()
