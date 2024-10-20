#!/bin/bash
# Job name:
#SBATCH --job-name=ipy
#
# Partition:
#SBATCH --partition=savio
#
# QoS:
#SBATCH --qos=nuclear_savio_normal
# #SBATCH --qos=savio_normal
# #SBATCH --qos=savio_debug
#
# Account:
#SBATCH --account=co_nuclear
# #SBATCH --account=fc_neutronics
#
# Processors:
#SBATCH --nodes=1
#
# Tasks per node (based on number of cores per node):
#SBATCH --ntasks-per-node=20

#
# Wall clock limit:
#SBATCH --time=24:00:00
##################################################################################################
module load python/2.7
export PYTHONIOENCODING=UTF-8
ipcluster start -n $SLURM_NTASKS &
sleep 45 # wait until all engines have successfully started
ipython ipy_test.py > job.pyout
