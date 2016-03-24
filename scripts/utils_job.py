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
import time
import math
import subprocess
import math
import sys
from collections import OrderedDict
import uuid

def get_python_job(job_id,
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
                   tbslas_output_dir=None):

    if not tbslas_output_dir:
        tbslas_output_dir = os.environ['TBSLAS_RESULT_DIR']
    if not os.path.exists(tbslas_output_dir):
        os.makedirs(tbslas_output_dir)

    TIMESTR       = time.strftime("%Y%m%d-%H%M%S")

    pyscript=\
"""
import os
from collections import OrderedDict
import utils
def run():
    ID=\'{PY_ID}\'
    prog=\'{PY_PROG}\'
    pn_list={PY_PN}
    tl_list={PY_TL}
    dp_list={PY_DP}
    cq_list={PY_CQ}
    ci_list={PY_CI}
    uf_list={PY_UF}
    np_list={PY_NP}
    nt_list={PY_NT}
    dt_list={PY_DT}
    tn_list={PY_TN}
    vs_list={PY_VS}
    mg_list={PY_MG}
    tt_list={PY_TT}

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
    utils.execute_commands(cmd_args, ID)

if __name__ == '__main__':
    tbslas_dir = os.environ['TBSLAS_RESULT_DIR']
    os.environ['TBSLAS_RESULT_DIR'] = \'{PY_OUTPUT_DIR}\'
    run()
    os.environ['TBSLAS_RESULT_DIR'] = tbslas_dir
"""
    output_dir='./'
    file_name= os.path.join(output_dir, job_id+'_'+TIMESTR+'_'+str(uuid.uuid4())+'.pyjob')
    file_handler = open(file_name,"w")
    global pyscript
    commands=pyscript.format(PY_ID=ID,
                             PY_PROG=prog,
                             PY_PN=pn_list,
                             PY_TL=tl_list,
                             PY_DP=dp_list,
                             PY_CQ=cq_list,
                             PY_CI=ci_list,
                             PY_UF=uf_list,
                             PY_NP=np_list,
                             PY_NT=nt_list,
                             PY_DT=dt_list,
                             PY_TN=tn_list,
                             PY_VS=vs_list,
                             PY_MG=mg_list,
                             PY_TT=tt_list,
                             PY_OUTPUT_DIR=tbslas_output_dir)

    file_handler.write(commands)
    file_handler.close()
    return file_name
