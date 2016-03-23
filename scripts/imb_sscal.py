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
    32,
    # 64
    ]
    num_threads = 10
    for nn in nn_list:
        total_time = '04:00:00'
        prog  = 'advection'
        # prog  = 'advdiff-ss'
        mg_list = [cnt      for cnt in range(1,4)]
        # mg_list = [3]
        num_pnts = 8**3
        num_steps = len(mg_list)
        np_list = [2*nn   for cnt in range(0,num_steps)]
        pn_list = [num_pnts        for cnt in range(0,num_steps)]
        tl_list = [1e-5            for cnt in range(0,num_steps)]
        dp_list = [10              for cnt in range(0,num_steps)]
        cq_list = [5               for cnt in range(0,num_steps)]
        ci_list = [False           for cnt in range(0,num_steps)]
        uf_list = [4               for cnt in range(0,num_steps)]
        nt_list = [num_threads for cnt in range(0,num_steps)]
        dt_list = [0.25            for cnt in range(0,num_steps)]
        tn_list = [1               for cnt in range(0,num_steps)]
        vs_list = [0               for cnt in range(0,num_steps)]
        tt_list = [6               for cnt in range(0,num_steps)]
        ID=prog+'-sscal-np-'+str(2*nn).zfill(5)
        
        print tbslas_res_dir
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
        sj.submit_job(job_name, nn, 2*nn, num_threads, total_time)

    os.environ['TBSLAS_RESULT_DIR'] = tbslas_dir

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
        sscal()





