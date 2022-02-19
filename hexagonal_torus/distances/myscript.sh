#!/bin/bash -l
#SBATCH -p high
#SBATCH -t 500:00:00

srun ./find_manifold_dim

