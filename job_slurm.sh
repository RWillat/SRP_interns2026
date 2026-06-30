#!/bin/bash

# Submit this script with: sbatch thefilename

#SBATCH --time=1:00:00 # walltime
#SBATCH --ntasks-per-node=1 # number of processor cores (i.e. tasks)
#SBATCH --nodes=2 # number of nodes
#SBATCH --mem=500m
#SBATCH --wckey edu_class # Project Code
#SBATCH -J "HPMR-GeomTest" # job name
#SBATCH --mail-user=rdwillat@umich.edu # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load use.moose  
# module load moose-apps        
# module load griffin-openmpi/2026.02.27-a7e65e5
module load griffin-openmpi

# mpirun -n 12 griffin-opt -i PNSource_test_RZ.i > run.log 2&1
mpiexec=/apps/local/openmpi/5.0.5-gcc13.2.0-container/bin/mpirun
griffinexec=/hpc-common/moose/containers/griffin-openmpi/2026.02.27-a7e65e5_ee1118e844308a213b63dfac1474f31ad5597d1e10d92be98e22aa8724c96338/bin/griffin-opt
#inputfile=GriffinTests/HPMR_OneSixth_Core_meshgenerator_tri.i
inputfile=GriffinTests/NEP_HPMR_meshgenerator_tri.i

$mpiexec -n 20 $griffinexec -i $inputfile # > run.log 2&1