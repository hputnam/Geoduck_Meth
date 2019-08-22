#!/bin/bash
## Job Name
#SBATCH --job-name=super-duck1
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes (We only get 1, so this is fixed)
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=10-100:00:00
## Memory per node
#SBATCH --mem=500G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr320@uw.edu
## Specify the working directory for this job
#SBATCH --workdir=/gscratch/srlab/sr320/analyses/0305





source /gscratch/srlab/programs/scripts/paths.sh


/gscratch/srlab/programs/supernova-2.0.0/supernova mkoutput \
        --asmdir=/gscratch/scrubbed/sr320/0220/Geoduck/outs/assembly \
        --outprefix=/gscratch/srlab/sr320/analyses/0305/duck4 \
        --style=pseudohap2