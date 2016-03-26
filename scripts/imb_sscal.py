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
import utils_job as job
import utils_submit_job as sj

def sscal():
    ############################################################################
    # STRONG SCALING
    ############################################################################
    tbslas_dir = os.environ['TBSLAS_RESULT_DIR']
    import time
    TIMESTR       = time.strftime("%Y%m%d-%H%M%S-")+str(time.time())
    tbslas_res_dir =  os.path.join(tbslas_dir,'imb-sscal-'+TIMESTR)
    if not os.path.exists(tbslas_res_dir):
        os.makedirs(tbslas_res_dir)

    nn_list = [
    1,
    2,
    4,
    8,
    16,
    # 32,
    # 64,
    # 128,
    # 256,
    # 512
    ]
    for num_nodes in nn_list:
        total_time = '00:30:00'
        num_threads = 7
        num_procs   = 4*num_nodes
        queue = 'micro'

        # prog  = 'advdiff-ss'
        merge_type_list = [1, 3, 2]
        for mg in merge_type_list:
            mg_list = [mg]
            prog  = 'advection'
            num_pnts = 8**5
            num_steps = len(mg_list)
            np_list = [num_procs       for cnt in range(0,num_steps)]
            nt_list = [num_threads for cnt in range(0,num_steps)]
            pn_list = [num_pnts        for cnt in range(0,num_steps)]
            tl_list = [1e-9            for cnt in range(0,num_steps)]
            dp_list = [9              for cnt in range(0,num_steps)]
            cq_list = [5               for cnt in range(0,num_steps)]
            ci_list = [False           for cnt in range(0,num_steps)]
            uf_list = [4               for cnt in range(0,num_steps)]
            dt_list = [0.25            for cnt in range(0,num_steps)]
            tn_list = [1               for cnt in range(0,num_steps)]
            vs_list = [0               for cnt in range(0,num_steps)]
            tt_list = [6               for cnt in range(0,num_steps)]
            ID=prog+'-sscal-np-'+str(num_procs).zfill(5)+'-mt-'+str(mg).zfill(2)

            job_name= job.get_python_job('imb-sscal',
                                         ID,
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
                                         tt_list,
                                         tbslas_res_dir)
            sj.submit_job(job_name, num_nodes, num_procs, num_threads, total_time, queue)

    os.environ['TBSLAS_RESULT_DIR'] = tbslas_dir

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
        sscal()
