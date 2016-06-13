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

# def generate_command_args(tl_list,\
#                           dt_list,\
#                           tn_list,\
#                           # de_list,\
#                           # q_list, \
#                           np_list,\
#                           nt_list,\
#                           use_cubic,\
#                           num_steps):

#     EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "advdiff-ss")
#     # EXEC = os.path.join(utils.TBSLAS_EXAMPLES_BIN_DIR, "advdiff-ss-tv")

#     # generate a dictionary data type of commands
#     cmd_args = OrderedDict()
#     cmd_id = 1;
#     for counter in range(0,num_steps):
#         ARGS    = ['-N'   , str(8**( math.ceil(math.log(np_list[counter],8))+1 ) ), \
#                    '-tol' , str(tl_list[counter]),                                  \
#                    '-q'   , str(14),\
#                    '-dt'  , str(dt_list[counter]),                                  \
#                    '-tn'  , str(tn_list[counter]),                                  \
#                    '-vsr', str(0),\
#                    '-omp' , str(nt_list[counter])]
#         if use_cubic:
#             ARGS = ARGS + ['-cubic', '1']
#         cmd_args[cmd_id] = utils.determine_command_prefix(np_list[counter]) + [EXEC] + ARGS
#         cmd_id = cmd_id + 1
#     return cmd_args

def conv_temporal():
    mpi_num_procs, omp_num_threads = utils.parse_args()
    num_steps = 8
    ############################################################################
    # TEST 1: TEMPORAL CONVERGENCE
    ############################################################################
    # prog          = 'advdiff-ss'
    # prog          = 'advdiff-ss-tv'
    prog          = 'advdiff-ss-tv-extrap'
    vtk_save_rate = 0
    mrg_type      = 3
    np            = mpi_num_procs
    num_pnts      = 8**(math.floor(math.log(np,8)+1))
    nt            = omp_num_threads

    # TREE TOLERANCE
    tl_factor = 1#0.1
    tl_init   = 1e-5
    tl_list   = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0,num_steps)]

    # TIME RESOLUTION
    dt_factor = 0.5
    # dt_init   = 0.1/16#0.5**5
    dt_init   = 0.5**0
    dt_list   = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0,num_steps)]

    # NUM TIME STEPS
    T_END     = 1.0
    tn_factor = 1.0/dt_factor
    tn_init   = T_END/dt_init
    tn_list   = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0,num_steps)]

    dp_list = [15            for cnt in range(0,num_steps)]
    cq_list = [14            for cnt in range(0,num_steps)]
    ci_list = [True          for cnt in range(0,num_steps)]
    uf_list = [2             for cnt in range(0,num_steps)]
    pn_list = [num_pnts      for cnt in range(0,num_steps)]
    np_list = [np            for cnt in range(0,num_steps)]
    nt_list = [nt            for cnt in range(0,num_steps)]
    mg_list = [mrg_type      for cnt in range(0,num_steps)]
    vs_list = [vtk_save_rate for cnt in range(0,num_steps)]
    tt_list = [11            for cnt in range(0,num_steps)]

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
    utils.execute_commands(cmd_args, 'temporal')

def conv_spatial():
    ############################################################################
    # TEST 2: SPATIAL CONVERGENCE
    ############################################################################
    use_cubic     = True
    # TREE TOLERANCE
    tl_factor = 0.1
    tl_init   = 1e-1
    tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0,num_steps)]

    # TIME RESOLUTION
    dt_factor = 1
    dt_init   = 1e-3
    dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0,num_steps)]

    # NUM TIME STEPS
    tn_factor = 1.0/dt_factor
    tn_init   = 1#T_END/dt_init
    tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0,num_steps)]

    # NUM MPI PROCESSES
    np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(tl_list,\
                                     dt_list,\
                                     tn_list,\
                                     # de_list,\
                                     # q_list, \
                                     np_list,\
                                     nt_list,\
                                     use_cubic,\
                                     num_steps)

    utils.execute_commands(cmd_args,'spatial')

def conv_temporal_spatial():
    ############################################################################
    # TEST 3: TEMPORAL/SPATIAL CONVERGENCE
    ############################################################################
    use_cubic     = True
    # TREE TOLERANCE
    tl_factor = 0.1
    tl_init   = 1e-1
    tl_list = [tl_init*math.pow(tl_factor,float(cnt)) for cnt in range(0,num_steps)]

    # TIME RESOLUTION
    dt_factor = 0.5
    dt_init   = 1
    dt_list = [dt_init*math.pow(dt_factor,float(cnt)) for cnt in range(0,num_steps)]

    # NUM TIME STEPS
    tn_factor = 1.0/dt_factor
    tn_init   = T_END/dt_init
    tn_list = [tn_init*math.pow(tn_factor,float(cnt)) for cnt in range(0,num_steps)]

    # NUM MPI PROCESSES
    np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    cmd_args = generate_command_args(tl_list,\
                                     dt_list,\
                                     tn_list,\
                                     # de_list,\
                                     # q_list, \
                                     np_list,\
                                     nt_list,\
                                     use_cubic,\
                                     num_steps)

    utils.execute_commands(cmd_args,'temporal-spatial')

def conv_temporal_spatial_long_time():
    ############################################################################
    # TEST 2: CONVERGENCE TEST FOR ADVECTION-DIFFUSION SOLVER
    ############################################################################
    mpi_num_procs, omp_num_threads = utils.parse_args()
    # prog          = 'advdiff-ss'
    # prog          = 'advdiff-ss-tv'
    prog          = 'advdiff-ss-tv-extrap'

    dt            = 0.0628
    vtk_save_rate = 0
    mrg_type      = 3
    np            = mpi_num_procs
    num_pnts      = 8**(math.floor(math.log(np,8)+1))
    nt            = omp_num_threads

    # UNIFORM
    # dp_list = [5    , 6    , 7    ]
    # cq_list = [3    , 3    , 3    ]
    # ci_list = [True , True , True ]
    # uf_list = [2    , 2    , 2    ]
    # dt_list = [dt   , dt/2 , dt/4 ]
    # tn_list = [100  , 200  , 400  ]
    # num_steps = len(dp_list)
    # tl_list = [1e-30         for cnt in range(0,num_steps)]

    # ADAPTIVE
    # LOW ORDER
    tl_list = [1e-02, 1e-03, 1e-04 ]
    dp_list = [15   , 15   , 15    ]
    cq_list = [3    , 3    , 3     ]
    ci_list = [True , True , True  ]
    uf_list = [2    , 2    , 2     ]
    dt_list = [dt   , dt/2 , dt/3  ]
    tn_list = [100  , 200  , 300   ]

    # HIGH ORDER
    # tl_list = [1e-03, 1e-04, 1e-05, 1e-6   ]
    # dp_list = [15   , 15   , 15   , 15     ]
    # cq_list = [14   , 14   , 14   , 14     ]
    # ci_list = [True , True , True , True   ]
    # uf_list = [2    , 2    , 2    , 2      ]
    # dt_list = [dt/4 , dt/8 , dt/16, dt/32 ]
    # tn_list = [400  , 800  , 1600 , 3200  ]


    num_steps = len(dp_list)
    pn_list = [num_pnts      for cnt in range(0,num_steps)]
    np_list = [np            for cnt in range(0,num_steps)]
    nt_list = [nt            for cnt in range(0,num_steps)]
    mg_list = [mrg_type      for cnt in range(0,num_steps)]
    vs_list = [vtk_save_rate for cnt in range(0,num_steps)]
    tt_list = [11            for cnt in range(0,num_steps)]

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
    utils.execute_commands(cmd_args, 'temporal-spatial-long-time')

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    tbslas_dir = os.environ['TBSLAS_RESULT_DIR']
    import time
    TIMESTR       = time.strftime("%Y%m%d-%H%M%S-")+str(time.time())
    os.environ['TBSLAS_RESULT_DIR'] =  os.path.join(tbslas_dir,'conv-advdiff-'+TIMESTR)
    if not os.path.exists(os.environ['TBSLAS_RESULT_DIR']):
        os.makedirs(os.environ['TBSLAS_RESULT_DIR'])

    # conv_temporal()
    # conv_spatial()
    # conv_temporal_spatial()
    conv_temporal_spatial_long_time()

    os.environ['TBSLAS_RESULT_DIR'] = tbslas_dir
