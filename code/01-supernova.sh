
#SBATCH --workdir=/gscratch/scrubbed/sr320/0220

source /gscratch/srlab/programs/scripts/paths.sh



/gscratch/srlab/programs/supernova-2.0.0/supernova run --id=Geoduck \
--fastqs=/gscratch/scrubbed/sr320/Chrom \
--sample=Geoduck-1,Geoduck-2,Geoduck-5,Geoduck-6



#SBATCH --workdir=/gscratch/srlab/sr320/analyses/0305


/gscratch/srlab/programs/supernova-2.0.0/supernova mkoutput \
        --asmdir=/gscratch/scrubbed/sr320/0220/Geoduck/outs/assembly \
        --outprefix=/gscratch/srlab/sr320/analyses/0305/duck4 \
        --style=pseudohap2
