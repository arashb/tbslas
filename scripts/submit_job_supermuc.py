import os

#template job script
commands=\
"""#! /usr/bin/ksh
#@ shell = /usr/bin/ksh
#@ job_type = MPICH
#@ initialdir={TBSLAS_RESULT_DIR}/scripts
#@ job_name = {JOB_ID}
#@ class = {QUEUE}
#@ node_usage = not_shared
#@ wall_clock_limit = {TOTAL_TIME}
#@ network.MPI = sn_all,,us,,
#@ notification = never
#@ output = {JOB_ID}.$(jobid).out
#@ error =  {JOB_ID}.$(jobid).err

#@ node = {NUM_NODES}
#@ total_tasks = {NUM_PROCS}
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

module load python/2.7_anaconda
module load git
module load fftw/serial/3.3
module load gcc/4.9

export TBSLAS_RESULT_DIR={TBSLAS_RESULT_DIR}
export TBSLAS_DIR={TBSLAS_DIR}

module unload mpi.ibm
module load mpi.intel
./.run_python.sh {JOB_id} {NUM_PROCS} {NUM_THREADS}\n"""

def get_job_file (num_nodes, num_procs, num_threads, queue, job_id, TIMESTR, total_time):

	print os.environ['HOME']
	file_handler = open(job_id[:-3]+'_'+TIMESTR+'.cmd',"w")

	global commands
	commands=commands.format(JOB_ID=job_id[:-3],\
	JOB_id=job_id,\
	QUEUE=queue,\
	TOTAL_TIME=total_time,\
	NUM_NODES=num_nodes,\
	NUM_PROCS=num_procs,\
	NUM_THREADS=num_threads,\
	TBSLAS_RESULT_DIR=os.environ['TBSLAS_RESULT_DIR'],\
	TBSLAS_DIR=os.environ['TBSLAS_DIR'])

	file_handler.write(commands)
	file_handler.close()

	return './'+job_id[:-3]+'_'+TIMESTR+'.cmd'

