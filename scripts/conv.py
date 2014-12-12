#!/usr/bin/python

#*************************************************************************
#Copyright (C) 2014 by Arash Bakhtiari
#You may not use this file except in compliance with the License.
#You obtain a copy of the License in the LICENSE file.

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.
#*************************************************************************
import os, subprocess, math

################################################################################
# CONSTRUCT THE EXECUTION COMMAND
################################################################################
PROGRAM  = "zalesak"
WORK_DIR = "examples/"
BIN_DIR  = os.path.join(WORK_DIR, "bin/")
RES_DIR  = os.path.join(WORK_DIR, "results/")

EXEC    = os.path.join(BIN_DIR, PROGRAM)


################################################################################
# MAIN
################################################################################
def main():
    dt = 0.1047 #math.pi/30
    toli = 1
    tolf = 9
    tol_list = [math.pow(0.1,x) for x in range(toli,tolf)]
    for tol in tol_list:
        ARGS    = ['-N', '8', '-tol', str(tol), '-dt', str(dt)]
        COMMAND = [EXEC] + ARGS
        print("BINARY: " +  str(COMMAND))
        p = subprocess.Popen(COMMAND, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            if line.startswith('TOL:'):
                print line

if __name__ == '__main__':
    main()
