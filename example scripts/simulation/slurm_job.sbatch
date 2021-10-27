#!/bin/bash
#SBATCH -J CPU_run
#SBATCH -c 24
#SBATCH --mail-type=END,FAIL
#SBATCH --output='/home/s1736159/MATLAB/slurm logs/slurm-%j.out'

module load mathworks/matlab_r2019a
module load nvidia/cuda-9.0

matlab -nodesktop -nosplash -nodisplay -r "batch_script_v2;"
