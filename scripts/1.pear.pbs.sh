#!/bin/bash
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N PEAR
#PBS -o PEAR.out
#PBS -e PEAR.err
#PBS -m n
#PBS -l nodes=1:thinnode:ppn=20
#PBS -l mem=100gb
#PBS -l walltime=20:00:00

set -e

module purge
module load tools
module load pear-academic/0.9.11


# Create output folder
out_dir="1_merged_reads"
mkdir --parents "$out_dir"
cd "$out_dir"


# Merge forward and reverse reads
pear \
--forward-fastq "S0_L001_R1_001.fastq.gz" \
--reverse-fastq "S0_L001_R2_001.fastq.gz" \
--output "LibName" \
--memory 100G \
--threads 20

# Delete unassembled and discarded reads
rm --force "LibName.discarded.fastq"
rm --force "LibName.unassembled.forward.fastq"
rm --force "LibName.unassembled.reverse.fastq"

echo "Paired reads merged successfully"


# Print debug info (optional, TORQUE only)
echo -e "\n\n\n---------- JOB INFO ----------\n\n"
qstat -f -1 "$PBS_JOBID"
