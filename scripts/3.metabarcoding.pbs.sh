#!/bin/bash
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N Metabarcoding
#PBS -o metabarcoding.out
#PBS -e metabarcoding.err
#PBS -m n
#PBS -l nodes=1:thinnode:ppn=20
#PBS -l mem=50gb
#PBS -l walltime=20:00:00

set -e

module purge
module load tools
module load anaconda2/4.4.0
module load anaconda3/2023.09-0
module load perl/5.36.1
module load ncbi-blast/2.15.0+


# Set paths
library='LibName'
assembled_fastq="1_merged_reads/$library.assembled.fastq"

blast_database="BLAST-North-Sea/12S/MIDORI2_UNIQ_NUC_GB261_srRNA"

taxonomy_dir="MIDORI2-GenBank261/QIIME/taxonomy"
taxonomy_file="$taxonomy_dir/expanded_ncbi_taxonomy.tsv"

blast_script_dir="scripts"
blast_assign_script="$blast_script_dir/taxonomy_assignment_BLAST.py"

# ------------------------------------------------------------------

# Demultiplex FASTQ file: remove tags
out_dir="2_demultiplex/$library"
mkdir --parents "$out_dir/trimmed_tags"
cd "$out_dir"

cutadapt \
--front "file:2_demultiplex/adapters/$library.tags.fasta" \
--overlap 7 \
--error-rate 0.1 \
--minimum-length 50 \
--maximum-length 350 \
--discard-untrimmed \
--rename '{id} sample={adapter_name}; tag_matches={match_sequence};' \
--cores 20 \
--output "$out_dir/trimmed_tags/$library.{name}.trimmed.fastq" \
"$assembled_fastq"

cat "$out_dir/trimmed_tags"/*.trimmed.fastq \
> "$out_dir/$library.no-tags.fastq"

echo "$library: tags trimmed"

# ------------------------------------------------------------------

# Demultiplex FASTQ file: remove primers
mkdir --parents "$out_dir/trimmed_primers"

cutadapt \
--front "file:2_demultiplex/adapters/$library.primers.fasta" \
--error-rate 0.1 \
--minimum-length 50 \
--maximum-length 350 \
--discard-untrimmed \
--rename '{header} primer={adapter_name}; primer_matches={match_sequence};' \
--cores 20 \
--output "$out_dir/trimmed_primers/$library.{name}.trimmed.fastq" \
"$out_dir/$library.no-tags.fastq"

cat "$out_dir/trimmed_primers"/*.trimmed.fastq \
> "$out_dir/$library.trimmed.fastq"

echo "$library: primers trimmed"

# ------------------------------------------------------------------

# Cluster identical reads
out_dir="3_obiuniq"
mkdir --parents "$out_dir"
cd "$out_dir"


obiuniq \
--merge sample \
--without-progress-bar \
"2_demultiplex/$library/$library.trimmed.fastq" \
> "$out_dir/$library.unique.fasta"

echo "$library: unique reads merged"


obistat \
--category-attribute count \
--without-progress-bar \
"$out_dir/$library.unique.fasta" \
| sort --numeric-sort --key 1 \
--output "$out_dir/$library.unique.stats.txt"

echo "$library: generated count statistics"

# ------------------------------------------------------------------

out_dir="4_obiclean"
mkdir --parents "$out_dir"
cd "$out_dir"

obigrep \
--without-progress-bar \
--lmin 20 \
--predicat 'count >= 2' \
"3_obiuniq/$library.unique.fasta" \
> "$out_dir/$library.pre-filtered.fasta"

echo "$library: pre-filtering completed"

# ------------------------------------------------------------------

obiclean \
--without-progress-bar \
--sample merged_sample \
--ratio 0.05 \
--head \
"$out_dir/$library.pre-filtered.fasta" \
> "$out_dir/$library.clean.fasta"

echo "$library: reads cleaned"

# ------------------------------------------------------------------

obiannotate \
--without-progress-bar \
--keep merged_sample \
--keep count \
--keep experiment \
--keep seq_length \
--keep forward_primer \
--keep reverse_primer \
"$out_dir/$library.clean.fasta" \
| obitab \
--without-progress-bar \
--no-definition \
> "$out_dir/$library.$target.tags.tsv"

gzip --best --force \
"$out_dir/$library.$target.tags.tsv"

rm --force "$out_dir/$library.$target.tags.tsv"

echo "$library: exported clusters' metadata"

# ------------------------------------------------------------------

# Run BLAST analysis
fasta_dir="4_obiclean"
fasta_file="$fasta_dir/$library.clean.fasta"

out_dir="4_BLAST/$library"
mkdir --parents "$out_dir"
cd "$out_dir"

blastn \
-evalue 1e-3 \
-perc_identity 97 \
-qcov_hsp_perc 97 \
-query "$fasta_file" \
-task megablast \
-db "$blast_database" \
-out "$out_dir/$library.raw.blast.tsv" \
-outfmt '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore staxids' \
-num_threads 20 \
-mt_mode 1

echo "$library: BLAST query executed"


# Refine BLAST results
python3 \
"$blast_assign_script" \
--hits_to_consider 5 \
--length_percentage 0.97 \
--percent_sway 0.5 \
--ncbi_nt \
--blast_database IGNORE \
--blast_file "$out_dir/$library.raw.blast.tsv" \
--output_dir "$out_dir" \
"$fasta_file" \
"$taxonomy_file"

echo "$library: OTU assigned to taxonomy"


# Print debug info (optional, TORQUE only)
echo -e "\n\n\n---------- JOB INFO ----------\n\n"
qstat -f -1 "$PBS_JOBID"
