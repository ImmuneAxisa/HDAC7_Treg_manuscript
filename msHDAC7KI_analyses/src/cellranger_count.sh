#!/bin/bash
#SBATCH --partition=general
#SBATCH -c 20 --nodes=1
#SBATCH --mem 110g


# path to fastq containing folders
FqPath=""
# Sample name
sample_name=""

MMT="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/refdata-gex-mm10-2020-A/"

cellranger count --id=$sample_name \
                        --transcriptome=$MMT \
                        --fastqs=$FqPath \
                        --sample=$sample_name > ${sample_name}.out 2>&1