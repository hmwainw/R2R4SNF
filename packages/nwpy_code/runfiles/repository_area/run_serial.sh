#!/bin/bash
# Job name:
#SBATCH --job-name=apc
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
# Tasks per node based on number of cores per node (example):
#SBATCH --ntasks-per-node=20
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=72:00:00
##################################################################################################
module load mpi4py
module load python/2.7.8 numpy scipy
python2 run.py