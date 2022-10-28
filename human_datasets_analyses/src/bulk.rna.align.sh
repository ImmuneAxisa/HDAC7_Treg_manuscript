#!/bin/bash

# Matthew R Lincoln rna.align.sh script
#
#
# As currently configured, this script runs on Ruddle. This can be altered by
# changing the paths defined in Section 1.
#
# The following steps are performed:
#   1. Check that all required components are present
#   2. QC of raw reads with FastQC
#   3. Adapter trimming with Trimmomatic
#   4. QC of trimmed reads with FastQC
#   5. Alignment of trimmed reads with STAR
#   6. Compilation of RNA-seq and alignment summary statistics by Picard tools
#   7. Filtering to remove optical duplicates
#   8. Compilation of summary statistics on filtered reads
#   9. Quantiate reads to gene level with htseq-count
#   10. Quantitate reads to gene level with RSEM


################################################################################
#############################   Section 1: Setup   #############################
################################################################################

# This section defines executable paths and locations of raw data. Any missing
# output directories are created.

missing=0 # Stores whether any critical components are missing

# project directory name
project_name= # project directory

# Raw data and working directories:
raw_data_dir="/home/pa326/project/raw_data/${project_name}"
base_dir="/home/pa326/project/${project_name}"
data_dir="${base_dir}/data"
scratch_dir="/home/pa326/scratch60/${project_name}"

# Output directories:
slurm_dir="${base_dir}/slurm"
logs_dir="${base_dir}/logs"
doc_dir="${base_dir}/doc"
results_dir="${base_dir}/results/rna"

# Scratch directories
pretrim_qc_dir="${scratch_dir}/rna/pretrim_fastqc"
posttrim_qc_dir="${scratch_dir}/rna/posttrim_fastqc"
trimmed_reads_dir="${scratch_dir}/rna/trimmed_reads"
alignment_dir="${scratch_dir}/rna/aligned_reads"
alignment_stats_dir="${scratch_dir}/rna/alignment_stats"
quantitation_dir="${scratch_dir}/rna/quantitation/star"
qc_stats_dir="${scratch_dir}/rna/qc_stats"

# Paths to program binaries:
fastqc_dir="/gpfs/ycga/home/pa326/bin/FastQC"
star_dir="/gpfs/ycga/home/pa326/bin/STAR-2.5.3a/bin/Linux_x86_64"
rsem_dir="/gpfs/ycga/home/pa326/bin/RSEM-1.3.0"

# Create any missing directories
directories=( $data_dir $logs_dir $doc_dir $slurm_dir $pretrim_qc_dir \
  $posttrim_qc_dir $trimmed_reads_dir $alignment_dir $alignment_stats_dir \
  $quantitation_dir/rsem $quantitation_dir/htseq $qc_stats_dir $results_dir)
for directory in ${directories[@]}; do
  if [ ! -d $directory ]; then
    echo "Creating directory ${directory}"
    mkdir -p $directory
  fi
done

# The sample table:
sample_table= # path to sample table, assuming tab separated fields
sample_table_header= # header for the table


# Adapter file and genome indices:
adapter_file="/ycga-gpfs/apps/hpc/Apps/Trimmomatic/0.33/Trimmomatic-0.33/adapters/NexteraPE-PE.fa"
star_index="/ycga-gpfs/project/ysm/hafler/mrl54/genomes/GRCh38.r88/STAR_index_primary_assembly/"
rsem_index="/ycga-gpfs/project/ysm/hafler/mrl54/genomes/GRCh38.r88/RSEM_index"
gene_annotation="/ycga-gpfs/project/ysm/hafler/mrl54/genomes/GRCh38.r88/Homo_sapiens.GRCh38.88.primary_assembly.gtf"

# Check that all essential components are present

essential_files=( $sample_table $adapter_file ${star_index}/SA \
                  ${rsem_index}.seq $gene_annotation ${fastqc_dir}/fastqc \
                  ${star_dir}/STAR ${rsem_dir}/rsem-calculate-expression )
for file in ${essential_files[@]}; do
  if [ ! -f $file ]; then
    echo "Cannot find ${file}"
    missing=1
  fi
done

# Die if any required components missing
if [ $missing = 1 ]; then
    echo "ERROR: Cannot find all required components."
    exit -1
fi


################################################################################
###################   Section 2: Align and Quantitate Reads   ##################
################################################################################

# This section submits multiple slurm scripts to perform quality control,
# alignment and quantitation

# Loop over samples and perform alignments:
while read $sample_table_header
do
  # Skip the header row and any samples that are not sequenced:
  [ $Sample_Number == "Sample_Number" ] && continue
  ([ $RNA_FastQ1 == "NA" ] || [ $RNA_FastQ2 == "NA" ]) && continue

  echo -n "Processing sample ${Sample_Name}..."

  # Make sure both input files exist and that analysis hasn't been done yet
  if [ ! -f ${raw_data_dir}/${RNA_FastQ1} ] || [ ! -f ${raw_data_dir}/${RNA_FastQ2} ]; then
    echo "could not find both input files; skipping."
    continue
  elif [ -f ${quantitation_dir}/rsem/${Sample_Name}.genes.results ]; then
    echo "results file found."
    continue
  else
    echo "creating Slurm script."
  fi

  # The script output file:
  output_script="${slurm_dir}/${Sample_Name}.rna.align.slurm"

  # Intermediate files:
  trim_pair1="${trimmed_reads_dir}/${Sample_Name}_pair_R1.fastq.gz"
  trim_unpair1="${trimmed_reads_dir}/${Sample_Name}_unpair_R1.fastq.gz"
  trim_pair2="${trimmed_reads_dir}/${Sample_Name}_pair_R2.fastq.gz"
  trim_unpair2="${trimmed_reads_dir}/${Sample_Name}_unpair_R2.fastq.gz"
  sorted_aligned_file="${alignment_dir}/${Sample_Name}_sortpos.bam"
  filtered_aligned_file="${alignment_dir}/${Sample_Name}_nodups.bam"
  transcript_aligned_file="${alignment_dir}/${Sample_Name}_transcriptome.bam"
  rnaseq_stats_file="${alignment_stats_dir}/${Sample_Name}_RNAseq_stats_genome.txt"
  alignment_stats_file="${alignment_stats_dir}/${Sample_Name}_alignment_stats_genome.txt"
  filtered_rnaseq_stats_file="${alignment_stats_dir}/${Sample_Name}_RNAseq_stats_genome_filtered.txt"
  filtered_alignment_stats_file="${alignment_stats_dir}/${Sample_Name}_alignment_stats_genome_filtered.txt"
  rnaseq_stats_transcript_file="${alignment_stats_dir}/${Sample_Name}_RNAseq_stats_transcriptome.txt"
  alignment_stats_transcript_file="${alignment_stats_dir}/${Sample_Name}_alignment_stats_transcriptome.txt"

  echo \
"#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=${Sample_Name}.align
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64gb
#SBATCH -o ${logs_dir}/${Sample_Name}.align.out
#SBATCH -e ${logs_dir}/${Sample_Name}.align.err

# Load required modules
module load SAMtools
module load GCC

export PATH=\$PATH:${fastqc_dir}:${star_dir}:${rsem_dir}

# Check MD5 sums: DISABLED
#status_1=\$(echo \"${RNA_MD5_1} ${raw_data_dir}/${RNA_FastQ1}\" | \\
 # md5sum -c | \\
  #cut -d ' ' -f 2)
#status_2=\$(echo \"${RNA_MD5_2} ${raw_data_dir}/${RNA_FastQ2}\" | \\
 # md5sum -c | \\
  #cut -d ' ' -f 2)

#if [ \$status_1 != "OK" ] || [ \$status_2 != "OK" ]; then
 # echo "ERROR: MD5 sums differ"
  #die
#fi

# Perform QC testing on raw reads with FastQC:
fastqc ${raw_data_dir}/${RNA_FastQ1} \\
       ${raw_data_dir}/${RNA_FastQ2} \\
       --outdir=$pretrim_qc_dir

# Trim reads for adapters and quality with Trimmomatic:
module load Trimmomatic
java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE \\
     -threads 40 \\
     ${raw_data_dir}/${RNA_FastQ1} \\
     ${raw_data_dir}/${RNA_FastQ2} \\
     $trim_pair1 \\
     $trim_unpair1 \\
     $trim_pair2 \\
     $trim_unpair2 \\
     ILLUMINACLIP:${adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:5

# Perform QC testing on trimmed reads with FastQC:
fastqc $trim_pair1 \\
       $trim_pair2 \\
       --outdir=$posttrim_qc_dir

# This is a bit convoluted:
# 1. Have STAR put its alignment, along with summary files, into subfolder of ${alignment_dir}/logs
# 2. Convert the SAM alignment into sorted BAM, place this in alignment directory
# 3. Delete original SAM file
# This has the effect of storing only one sorted BAM alignment, and keeps the logs in a consistent place
# Do alignment with STAR:
if [ ! -d ${alignment_dir}/logs/${Sample_Name} ]; then
  mkdir -p ${alignment_dir}/logs/${Sample_Name}
fi

# with two-pass mapping per sample (consider using multi-sample two-pass) and quantitation to transcriptome (for RSEM):
STAR --runThreadN 40 \\
     --genomeDir $star_index \\
     --readFilesIn $trim_pair1 $trim_pair2 \\
     --readFilesCommand zcat \\
     --outFileNamePrefix ${alignment_dir}/logs/${Sample_Name}/ \\
     --outSAMtype BAM SortedByCoordinate \\
     --outSAMstrandField intronMotif \\
     --outFilterIntronMotifs RemoveNoncanonical \\
     --outFilterType BySJout \\
     --outFilterMultimapNmax 20 \\
     --alignSJoverhangMin 8 \\
     --alignSJDBoverhangMin 1 \\
     --outFilterMismatchNmax 999 \\
     --outFilterMismatchNoverReadLmax 0.04 \\
     --alignIntronMin 20 \\
     --alignIntronMax 1000000 \\
     --alignMatesGapMax 1000000 \\
     --quantMode TranscriptomeSAM \\
     --twopassMode Basic

if grep 'ALL DONE!' ${alignment_dir}/logs/${Sample_Name}/Log.progress.out; then
  mv ${alignment_dir}/logs/${Sample_Name}/Log.final.out ${alignment_dir}/logs/${Sample_Name}.star.summary.log
  mv ${alignment_dir}/logs/${Sample_Name}/Aligned.sortedByCoord.out.bam \\
    $sorted_aligned_file
  mv ${alignment_dir}/logs/${Sample_Name}/Aligned.toTranscriptome.out.bam \\
    $transcript_aligned_file
  rm -rf ${alignment_dir}/logs/${Sample_Name}
fi

module load picard

# Collect alignment statistics for reads aligned to the genome:
java -jar \$EBROOTPICARD/picard.jar CollectRnaSeqMetrics \\
     I=$sorted_aligned_file \\
     O=$rnaseq_stats_file \\
     REF_FLAT=/ycga-gpfs/project/ysm/hafler/mrl54/genomes/GRCh38.r88/Homo_sapiens.GRCh38.88.refFlat \\
     STRAND=NONE
java -jar \$EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \\
     R=/ycga-gpfs/project/ysm/hafler/mrl54/genomes/GRCh38.r88/Homo_sapiens.GRCh38.dna.primary_assembly.fa \\
     I=$sorted_aligned_file \\
     O=$alignment_stats_file

# Filter optical duplicates:
java -jar \$EBROOTPICARD/picard.jar MarkDuplicates \\
     I=$sorted_aligned_file \\
     O=$filtered_aligned_file \\
     M=${logs_dir}/${Sample_Name}_star_duplicates.log \\
     REMOVE_DUPLICATES=TRUE \\
     VALIDATION_STRINGENCY=SILENT
# rm $sorted_aligned_file

# Collect post-filtering quality statistics for reads aligned to the genome:
java -jar \$EBROOTPICARD/picard.jar CollectRnaSeqMetrics \\
     I=$filtered_aligned_file \\
     O=$filtered_rnaseq_stats_file \\
     REF_FLAT=/ycga-gpfs/project/ysm/hafler/mrl54/genomes/GRCh38.r88/Homo_sapiens.GRCh38.88.refFlat \\
     STRAND=NONE
java -jar \$EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \\
     R=/ycga-gpfs/project/ysm/hafler/mrl54/genomes/GRCh38.r88/Homo_sapiens.GRCh38.dna.primary_assembly.fa \\
     I=$filtered_aligned_file \\
     O=$filtered_alignment_stats_file

# Quantitate reads to gene level with htseq-count:
module load GCC
htseq-count -f bam -r pos -s no \\
  $filtered_aligned_file \\
  $gene_annotation > \\
  ${quantitation_dir}/htseq/${Sample_Name}.gtf

# Collect alignment statistics for reads aligned to the transcriptome:
# java -jar \$EBROOTPICARD/picard.jar CollectRnaSeqMetrics \\
#      I=$transcript_aligned_file \\
#      O=$rnaseq_stats_transcript_file \\
#      REF_FLAT=/ycga-gpfs/project/ysm/hafler/mrl54/genomes/GRCh38.r88/Homo_sapiens.GRCh38.88.refFlat \\
#      STRAND=NONE
# java -jar \$EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \\
#      R=/ycga-gpfs/project/ysm/hafler/mrl54/genomes/GRCh38.r88/Homo_sapiens.GRCh38.dna.primary_assembly.fa \\
#      I=$transcript_aligned_file \\
#      O=$alignment_stats_transcript_file

# Quantitate STAR transcriptome-aligned reads to genes with RSEM:
module load RSEM
rsem-calculate-expression \\
  --no-bam-output \\
  -p 40 \\
  --bam \\
  --paired-end \\
  $transcript_aligned_file \\
  $rsem_index \\
  ${quantitation_dir}/rsem/${Sample_Name}" > $output_script

  sbatch $output_script
done < $sample_table

# Alternatively, can use UCSC refFlat file, converted to Ensembl chromosome names:
# cat /home/mrl54/project/genomes/GRCh38.r88/hg38.refFlat.txt | \
#   awk 'BEGIN { OFS="\t" } { sub(/chr/, "", $3); print }' > \
#   hg38.ensembl.chrs.refFlat.txt

# Subsequent steps in this script require jobs above to be completed:
exit

################################################################################
################   Section 3: Extract Sequencing QC statistics   ###############
################################################################################

mkdir -p ${results_dir}/qc_stats/

# Extract raw qc data from the fastqc archives:
while read $sample_table_header
do
  # Skip the header row and any samples that are not sequenced:
  [ $Sample_Number == "Sample_Number" ] && continue
  ([ $RNA_FastQ1 == "NA" ] || [ $RNA_FastQ2 == "NA" ]) && continue

  pretrim_r1_zip=`echo ${pretrim_qc_dir}/${RNA_FastQ1} | sed 's/\.fastq\.gz$/_fastqc.zip/'`
  pretrim_r2_zip=`echo ${pretrim_qc_dir}/${RNA_FastQ2} | sed 's/\.fastq\.gz$/_fastqc.zip/'`
  posttrim_r1_zip=`echo ${posttrim_qc_dir}/${RNA_FastQ2} | sed 's/\.fastq\.gz$/_fastqc.zip/'`
  posttrim_r2_zip=`echo ${posttrim_qc_dir}/${RNA_FastQ2} | sed 's/\.fastq\.gz$/_fastqc.zip/'`

  unzip -p $pretrim_r1_zip */fastqc_data.txt > \
    ${qc_stats_dir}/${Sample_Name}_R1.pretrim.fastqc.data.txt
  unzip -p $pretrim_r2_zip */fastqc_data.txt > \
    ${qc_stats_dir}/${Sample_Name}_R2.pretrim.fastqc.data.txt
  unzip -p ${posttrim_qc_dir}/${Sample_Name}_pair_R1_fastqc.zip */fastqc_data.txt > \
    ${qc_stats_dir}/${Sample_Name}_R1.posttrim.fastqc.data.txt
  unzip -p ${posttrim_qc_dir}/${Sample_Name}_pair_R2_fastqc.zip */fastqc_data.txt > \
    ${qc_stats_dir}/${Sample_Name}_R2.posttrim.fastqc.data.txt
done < $sample_table

# Sequence quality per base:
seq_qual_base=${results_dir}/qc_stats/sequence.quality.per.base.txt
echo -e "Sample\tPhase\tRead\tBase\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile" > \
  $seq_qual_base
while read $sample_table_header
do
  # Skip the header row and any samples that are not sequenced:
  [ $Sample_Number == "Sample_Number" ] && continue
  ([ $RNA_FastQ1 == "NA" ] || [ $RNA_FastQ2 == "NA" ]) && continue

  cat ${qc_stats_dir}/${Sample_Name}_R1.pretrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Per base sequence quality/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R1",$0 }' >> \
      $seq_qual_base
  cat ${qc_stats_dir}/${Sample_Name}_R2.pretrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Per base sequence quality/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R2",$0 }' >> \
      $seq_qual_base
  cat ${qc_stats_dir}/${Sample_Name}_R1.posttrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Per base sequence quality/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R1",$0 }' >> \
      $seq_qual_base
  cat ${qc_stats_dir}/${Sample_Name}_R2.posttrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Per base sequence quality/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R2",$0 }' >> \
      $seq_qual_base
done < $sample_table

# Sequence duplication levels:
seq_dup_level=${results_dir}/qc_stats/sequence.duplication.level.txt
echo -e "Sample\tPhase\tRead\tDuplication Level\tPercentage of deduplicated\tPercentage of total" > \
  $seq_dup_level
while read $sample_table_header
do
  # Skip the header row and any samples that are not sequenced:
  [ $Sample_Number == "Sample_Number" ] && continue
  ([ $RNA_FastQ1 == "NA" ] || [ $RNA_FastQ2 == "NA" ]) && continue

  cat ${qc_stats_dir}/${Sample_Name}_R1.pretrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Sequence Duplication Levels/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R1",$0 }' >> \
      $seq_dup_level
  cat ${qc_stats_dir}/${Sample_Name}_R2.pretrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Sequence Duplication Levels/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R2",$0 }' >> \
      $seq_dup_level
  cat ${qc_stats_dir}/${Sample_Name}_R1.posttrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Sequence Duplication Levels/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R1",$0 }' >> \
      $seq_dup_level
  cat ${qc_stats_dir}/${Sample_Name}_R2.posttrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Sequence Duplication Levels/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R2",$0 }' >> \
      $seq_dup_level
done < $sample_table

# GC content levels:
gc_content_level=${results_dir}/qc_stats/gc.content.level.txt
echo -e "Sample\tPhase\tRead\tGC Content\tCount" > \
  $gc_content_level
while read $sample_table_header
do
  # Skip the header row and any samples that are not sequenced:
  [ $Sample_Number == "Sample_Number" ] && continue
  ([ $RNA_FastQ1 == "NA" ] || [ $RNA_FastQ2 == "NA" ]) && continue

  cat ${qc_stats_dir}/${Sample_Name}_R1.pretrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Per sequence GC content/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R1",$0 }' >> \
      $gc_content_level
  cat ${qc_stats_dir}/${Sample_Name}_R2.pretrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Per sequence GC content/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R2",$0 }' >> \
      $gc_content_level
  cat ${qc_stats_dir}/${Sample_Name}_R1.posttrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Per sequence GC content/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R1",$0 }' >> \
      $gc_content_level
  cat ${qc_stats_dir}/${Sample_Name}_R2.posttrim.fastqc.data.txt | \
    awk -v sample=$Sample_Name 'BEGIN { OFS="\t" }
      /^>>Per sequence GC content/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R2",$0 }' >> \
      $gc_content_level
done < $sample_table


################################################################################
######################   Section 4: Alignment statistics   #####################
################################################################################

alignment_stats=${results_dir}/qc_stats/star.alignment.stats.txt
echo -e "Sample\tInput\tUnique\tMultiple\tToo many" > \
  $alignment_stats
while read $sample_table_header
do
  # Skip the header row and any samples that are not sequenced:
  [ $Sample_Number == "Sample_Number" ] && continue
  ([ $RNA_FastQ1 == "NA" ] || [ $RNA_FastQ2 == "NA" ]) && continue

  cat ${alignment_dir}/logs/${Sample_Name}.star.summary.log | awk -v sample=$Sample_Name \
    'BEGIN{ OFS="\t";
            input=0;
            unique=0;
            multiple=0;
            too_many=0; }
     /Number of input reads/ { input=$6 }
     /Uniquely mapped reads number/ { unique=$6 }
     /Number of reads mapped to multiple loci/ { multiple=$9 }
     /Number of reads mapped to too many loci/ { too_many=$10 }
     END{ print sample,input,unique,multiple,too_many }' >> \
    $alignment_stats
done < $sample_table


################################################################################
###################   Section 5: Post-alignment statistics   ###################
################################################################################

# Collect picard alignment statistics
picard_alignment=${results_dir}/qc_stats/picard.alignment.stats.txt
echo -e "Sample\tPhase\tClass\tCount" > \
  $picard_alignment
while read $sample_table_header
do
  # Skip the header row and any samples that are not sequenced:
  [ $Sample_Number == "Sample_Number" ] && continue
  ([ $RNA_FastQ1 == "NA" ] || [ $RNA_FastQ2 == "NA" ]) && continue

  cat ${alignment_stats_dir}/${Sample_Name}_alignment_stats_genome.txt | \
  awk -v sample=${Sample_Name} \
    'BEGIN{ OFS="\t" }
     /^PAIR/ { print sample,"Pre-filter","Total Reads",$2
               print sample,"Pre-filter","PF Reads",$3
               print sample,"Pre-filter","PF Reads Aligned",$6
               print sample,"Pre-filter","PF HQ Aligned Reads",$9
               print sample,"Pre-filter","Reads Aligned in Pairs",$17
             }' >> $picard_alignment

  cat ${alignment_stats_dir}/${Sample_Name}_alignment_stats_genome_filtered.txt | \
  awk -v sample=${Sample_Name} \
    'BEGIN{ OFS="\t" }
     /^PAIR/ { print sample,"Post-filter","Total Reads",$2
               print sample,"Post-filter","PF Reads",$3
               print sample,"Post-filter","PF Reads Aligned",$6
               print sample,"Post-filter","PF HQ Aligned Reads",$9
               print sample,"Post-filter","Reads Aligned in Pairs",$17
             }' >> $picard_alignment
done < $sample_table

# Collect picard RNA-seq statistics
picard_rnaseq=${results_dir}/qc_stats/picard.rnaseq.stats.txt
echo -e "Sample\tPhase\tClass\tCount" > \
  $picard_rnaseq
while read $sample_table_header
do
  # Skip the header row and any samples that are not sequenced:
  [ $Sample_Number == "Sample_Number" ] && continue
  ([ $RNA_FastQ1 == "NA" ] || [ $RNA_FastQ2 == "NA" ]) && continue

  cat ${alignment_stats_dir}/${Sample_Name}_RNAseq_stats_genome.txt |
    grep -A 1 ^PF_BASES |
    grep -v ^PF_BASES |
    awk -v sample=${Sample_Name} \
      'BEGIN { OFS="\t"}
             { print sample,"Pre-filter","PF Bases",$1
               print sample,"Pre-filter","PF Aligned Bases",$2
               print sample,"Pre-filter","Ribosomal Bases",$3
               print sample,"Pre-filter","Coding Bases",$4
               print sample,"Pre-filter","UTR Bases",$5
               print sample,"Pre-filter","Intronic Bases",$6
               print sample,"Pre-filter","Intergenic Bases",$7
             }' >> $picard_rnaseq

  cat ${alignment_stats_dir}/${Sample_Name}_RNAseq_stats_genome_filtered.txt |
    grep -A 1 ^PF_BASES |
    grep -v ^PF_BASES |
    awk -v sample=${Sample_Name} \
      'BEGIN { OFS="\t"}
             { print sample,"Post-filter","PF Bases",$1
               print sample,"Post-filter","PF Aligned Bases",$2
               print sample,"Post-filter","Ribosomal Bases",$3
               print sample,"Post-filter","Coding Bases",$4
               print sample,"Post-filter","UTR Bases",$5
               print sample,"Post-filter","Intronic Bases",$6
               print sample,"Post-filter","Intergenic Bases",$7
             }' >> $picard_rnaseq
done < $sample_table

################################################################################
##########################   Section 6: Count Matrix   #########################
################################################################################

# We will construct the count matrices in R, but to do this, we need to copy the
# quantitation files to a location that can be read in R

mkdir -p ${results_dir}/quantitation/rsem \
         ${results_dir}/quantitation/htseq-count

while read $sample_table_header
do
  # Skip the header row and any samples that are not sequenced:
  [ $Sample_Number == "Sample_Number" ] && continue
  ([ $RNA_FastQ1 == "NA" ] || [ $RNA_FastQ2 == "NA" ]) && continue

  # Copy RSEM results to results folder; these can be imported into R with
  # tximport()
  cp ${quantitation_dir}/rsem/${Sample_Name}.genes.results \
    ${results_dir}/quantitation/rsem

  # Copy htseq-count results to results folder; these can be imported into R
  # with DESeqDataSetFromHTSeqCount()
  cp ${quantitation_dir}/htseq/${Sample_Name}.gtf \
    ${results_dir}/quantitation/htseq-count
done < $sample_table
exit
