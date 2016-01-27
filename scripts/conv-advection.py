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


def conv_temporal():
    ############################################################################
    # TEMPORAL CONVERGENCE TEST FOR ADVECTION
    ############################################################################
    mpi_num_procs, omp_num_threads = utils.parse_args()
    prog = 'advection'
    num_steps = 8

    ##############################
    # TREE TOLERANCE
    ##############################
    tl_fact = 1
    tl_init = 1e-9
    tl_list = [tl_init*math.pow(tl_fact,float(cnt)) for cnt in range(0,num_steps)]

    ##############################
    # TIME RESOLUTION
    ##############################
    dt_fact = 0.5
    dt_init = 1
    dt_list = [dt_init*math.pow(dt_fact,float(cnt)) for cnt in range(0,num_steps)]

    T_END   = 1.0
    tn_fact = 1.0/dt_fact
    tn_init = T_END/dt_init
    tn_list = [tn_init*math.pow(tn_fact,float(cnt)) for cnt in range(0,num_steps)]

    ##############################
    # TREE DEPTH/POINTS
    ##############################
    dp_list  = [15 for cnt in range(0,num_steps)]

    num_pnts = 8**(math.floor(math.log(mpi_num_procs,8)+1))
    pn_list  = [num_pnts      for cnt in range(0,num_steps)]

    ##############################
    # PARALLEL
    ##############################
    mpi_num_procs = mpi_num_procs
    np_list = [mpi_num_procs for cnt in range(0,num_steps)]

    nt = omp_num_threads
    nt_list = [nt for cnt in range(0,num_steps)]

    mrg_type = 3
    mg_list = [mrg_type for cnt in range(0,num_steps)]

    ##############################
    # CHEBYSHEV/CUBIC INTERPOLATION
    ##############################
    cq_list = [14 for cnt in range(0, num_steps)]
    ci_list = [True for cnt in range(0, num_steps)]
    uf_list = [4     for cnt in range(0, num_steps)]

    ##############################
    # VISUALIZATION
    ##############################
    vtk_save_rate = 10
    vs_list = [vtk_save_rate for cnt in range(0,num_steps)]

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
    utils.execute_commands(cmd_args, 'temporal')

def conv_spatial():
    ############################################################################
    # TEMPORAL CONVERGENCE TEST FOR ADVECTION
    ############################################################################
    mpi_num_procs, omp_num_threads = utils.parse_args()
    prog = 'advection'
    num_steps = 8

    ##############################
    # TREE TOLERANCE
    ##############################
    tl_fact = 0.1
    tl_init = 1e-1
    tl_list = [tl_init*math.pow(tl_fact,float(cnt)) for cnt in range(0,num_steps)]

    ##############################
    # TIME RESOLUTION
    ##############################
    dt_fact = 1
    dt_init = 1e-3
    dt_list = [dt_init*math.pow(dt_fact,float(cnt)) for cnt in range(0,num_steps)]

    T_END   = 1.0
    tn_fact = 1.0/dt_fact
    tn_init = 1#T_END/dt_init
    tn_list = [tn_init*math.pow(tn_fact,float(cnt)) for cnt in range(0,num_steps)]

    ##############################
    # TREE DEPTH/POINTS
    ##############################
    dp_list  = [15 for cnt in range(0,num_steps)]

    num_pnts = 8**(math.floor(math.log(mpi_num_procs,8)+1))
    pn_list  = [num_pnts      for cnt in range(0,num_steps)]

    ##############################
    # PARALLEL
    ##############################
    mpi_num_procs = mpi_num_procs
    np_list = [mpi_num_procs for cnt in range(0,num_steps)]

    nt = omp_num_threads
    nt_list = [nt for cnt in range(0,num_steps)]

    mrg_type = 3
    mg_list = [mrg_type for cnt in range(0,num_steps)]

    ##############################
    # CHEBYSHEV/CUBIC INTERPOLATION
    ##############################
    cq_list = [14 for cnt in range(0, num_steps)]
    ci_list = [True for cnt in range(0, num_steps)]
    uf_list = [4     for cnt in range(0, num_steps)]

    ##############################
    # VISUALIZATION
    ##############################
    vtk_save_rate = 10
    vs_list = [vtk_save_rate for cnt in range(0,num_steps)]

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
    utils.execute_commands(cmd_args, 'temporal')

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    ############################################################################
    # TEST 1: TEMPORAL ERROR
    ############################################################################
    conv_temporal()

    ############################################################################
    # TEST 2: SPATIAL ERROR
    ############################################################################
    # conv_spatial()

    # ############################################################################
    # # TEST 3: TEMPORAL/SPATIAL ERROR
    # ############################################################################
    # tl_fact = 0.1
    # tl_init = 1e-1
    # tl_list = [tl_init*math.pow(tl_fact,float(cnt)) for cnt in range(0,num_steps)]

    # dt_fact = 0.5
    # dt_init = 1
    # dt_list = [dt_init*math.pow(dt_fact,float(cnt)) for cnt in range(0,num_steps)]

    # tn_fact = 1.0/dt_fact
    # tn_init = T_END/dt_init
    # tn_list = [tn_init*math.pow(tn_fact,float(cnt)) for cnt in range(0,num_steps)]

    # # NUM MPI PROCESSES
    # np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # # NUM OMP THREADS
    # nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    # cmd_args = generate_command_args(tl_list,\
    #                                  dt_list,\
    #                                  tn_list,\
    #                                  # de_list,\
    #                                  # q_list, \
    #                                  np_list,\
    #                                  nt_list,\
    #                                  num_steps)

    # utils.execute_commands(cmd_args, 'table3')
