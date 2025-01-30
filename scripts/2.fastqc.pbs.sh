#!/bin/bash
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N FastQC
#PBS -o FastQC.out
#PBS -e FastQC.err
#PBS -m n
#PBS -l nodes=1:thinnode:ppn=10
#PBS -l mem=50gb
#PBS -l walltime=5:00:00

set -e

module purge
module load tools
module load perl/5.36.1
module load openjdk/22
module load fastqc/0.11.9


# Set input and output folders
fastq_dir="1_merged_reads"
fastq_dir=$(realpath "$fastq_dir")

out_dir="1_QC"
out_dir=$(realpath "$out_dir")

mkdir --parents "$out_dir"
cd "$out_dir"


# Check quality of assembled reads
fastqc \
--outdir "$out_dir" \
--threads 10 \
--quiet \
--dir "$out_dir" \
"$fastq_dir/LibName.assembled.fastq"

echo "FastQC report generated successfully"


# Print debug info (optional, TORQUE only)
echo -e "\n\n\n---------- JOB INFO ----------\n\n"
qstat -f -1 "$PBS_JOBID"
