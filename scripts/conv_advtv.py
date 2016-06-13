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
    # prog = 'field-set'
    prog = 'advtvextrap'
    num_steps = 10

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
    vtk_save_rate = 0
    vs_list = [vtk_save_rate for cnt in range(0,num_steps)]

    test = 10
    tt_list = [test for cnt in range(0,num_steps)]

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
    # TEMPORAL CONVERGENCE TEST FOR ADVECTION
    ############################################################################
    mpi_num_procs, omp_num_threads = utils.parse_args()
    prog = 'advtv'
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
    cq_list = [14   for cnt in range(0, num_steps)]
    ci_list = [True for cnt in range(0, num_steps)]
    uf_list = [4    for cnt in range(0, num_steps)]

    ##############################
    # VISUALIZATION
    ##############################
    vtk_save_rate = 0
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
    utils.execute_commands(cmd_args, 'spatial')

def conv_temporal_spatial():
    mpi_num_procs, omp_num_threads = utils.parse_args()
    prog = 'advtv'
    num_steps = 8

    ############################################################################
    # TEST 3: TEMPORAL/SPATIAL ERROR
    ############################################################################
    tl_fact = 0.1
    tl_init = 1e-1
    tl_list = [tl_init*math.pow(tl_fact,float(cnt)) for cnt in range(0,num_steps)]

    dt_fact = 0.5
    dt_init = 1
    dt_list = [dt_init*math.pow(dt_fact,float(cnt)) for cnt in range(0,num_steps)]

    tn_fact = 1.0/dt_fact
    tn_init = T_END/dt_init
    tn_list = [tn_init*math.pow(tn_fact,float(cnt)) for cnt in range(0,num_steps)]

    # NUM MPI PROCESSES
    np_list = [mpi_num_procs  for cnt in range(0, num_steps)]

    # NUM OMP THREADS
    nt_list = [omp_num_threads for cnt in range(0, num_steps)]

    # cmd_args = generate_command_args(tl_list,\
    #                                  dt_list,\
    #                                  tn_list,\
    #                                  # de_list,\
    #                                  # q_list, \
    #                                  np_list,\
    #                                  nt_list,\
    #                                  num_steps)

    # utils.execute_commands(cmd_args, 'temporal-spatial')

def conv_temporal_spatial_long_time():
    ############################################################################
    # TEST 2: CONVERGENCE TEST FOR ADVECTION
    ############################################################################
    mpi_num_procs, omp_num_threads = utils.parse_args()
    prog          = 'advtv'
    prog          = 'advtvextrap'
    dt            = 0.0628
    vsr           = 0
    mrg_type      = 3
    np            = mpi_num_procs
    num_pnts      = 8**(math.floor(math.log(np,8)+1))
    nt            = omp_num_threads

    # UNIFORM
    # dp_list = [5    , 6    , 7    ]#, 5     , 6     , 7    ]
    # cq_list = [3    , 3    , 3    ]#, 3     , 3     , 3    ]
    # ci_list = [False , False , False ]#, False , False , False]
    # uf_list = [2    , 2    , 2    ]#, 2     , 2     , 2    ]
    # dt_list = [dt   , dt/2 , dt/4 ]#, dt    , dt/2  , dt/4 ]
    # tn_list = [100  , 200  , 400  ]#, 100   , 200   , 400  ]
    # num_steps = len(dp_list)
    # tl_list = [1e-30         for cnt in range(0,num_steps)]

    # ADAPTIVE
    # LOW ORDER
    tl_list = [1e-02, 1e-03, 1e-04 ]
    dp_list = [15   , 15   , 15    ]
    cq_list = [3    , 3    , 3     ]
    ci_list = [True , True , True  ]
    uf_list = [2    , 2    , 2     ]
    dt_list = [dt   , dt/2 , dt/4  ]
    tn_list = [100  , 200  , 400   ]

    # HIGH ORDER
    # tl_list = [1e-05, 1e-06, 1e-07, 1e-8]#, 1e-02 , 1e-03 , 1e-04]
    # dp_list = [15   , 15   , 15   , 15  ]#, 15    , 15    , 15   ]
    # cq_list = [14   , 14   , 14   , 14  ]#, 3     , 3     , 3    ]
    # ci_list = [True , True , True , True]#, False , False , False]
    # uf_list = [2    , 2    , 2    , 2   ]#, 2     , 2     , 2    ]
    # dt_list = [dt/4 , dt/8 , dt/16 , dt/32]#, dt    , dt/2  , dt/4 ]
    # tn_list = [400  , 800  , 1600  , 3200]#, 100   , 200   , 400  ]

    num_steps = len(dp_list)
    pn_list = [num_pnts      for cnt in range(0,num_steps)]
    np_list = [np            for cnt in range(0,num_steps)]
    nt_list = [nt            for cnt in range(0,num_steps)]
    mg_list = [mrg_type      for cnt in range(0,num_steps)]
    vs_list = [vsr           for cnt in range(0,num_steps)]
    tt_list = [10            for cnt in range(0,num_steps)]

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
    os.environ['TBSLAS_RESULT_DIR'] =  os.path.join(tbslas_dir,'conv-advection-'+TIMESTR)
    if not os.path.exists(os.environ['TBSLAS_RESULT_DIR']):
        os.makedirs(os.environ['TBSLAS_RESULT_DIR'])

    # conv_temporal()
    # conv_spatial()
    # conv_temporal_spatial()
    conv_temporal_spatial_long_time()

    os.environ['TBSLAS_RESULT_DIR'] = tbslas_dirl
