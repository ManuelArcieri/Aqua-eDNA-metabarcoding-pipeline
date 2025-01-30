#!/bin/bash
#PBS -W group_list=dp_00007 -A dp_00007
#PBS -N Make_BLAST_database
#PBS -o make-blast-database.out
#PBS -e make-blast-database.err
#PBS -m n
#PBS -l nodes=1:thinnode:ppn=1
#PBS -l mem=50gb
#PBS -l walltime=1:00:00

set -e

module purge
module load tools
module load perl/5.36.1
module load python/3.12.0
module load ncbi-blast/2.16.0+
module load seqkit/2.4.0


# Directory and file prefix of the MIDORI database (must be in QIIME format)
midori_dir="sample_data"
midori_dir=$(realpath "$midori_dir")
prefix="MIDORI2_UNIQ_NUC_GB263_srRNA_QIIME"

# Table containing a list of acceptable taxa IDs
acceptable_species_file="acceptable_species.tsv"
acceptable_species_file=$(realpath "$acceptable_species_file")

# Output directory for taxonomic information
taxonomy_dir="taxonomy"
taxonomy_dir=$(realpath "$taxonomy_dir")

# Folder containing "genbank_nodes_and_names_to_taxonomy.py"
blast_script_dir="scripts"
blast_script_dir=$(realpath "$blast_script_dir")
blast_taxonomy_script="$blast_script_dir/genbank_nodes_and_names_to_taxonomy.py"

# Output directory for the BLAST database
blast_out_dir="blast"
mkdir --parents "$blast_out_dir"
blast_out_dir=$(realpath "$blast_out_dir")


# Write a list of all haplotypes and their respective taxa IDs
cd "$midori_dir"

zcat "$prefix.taxon.gz" \
| awk \
--field-separator ';' \
'{print $1, $NF}' \
| awk \
--field-separator ' ' \
'{print $1, $NF}' \
| awk \
'{
    split($1, a, ".")
    first_field = a[1] "." a[2]

    split($2, b, "_")
    second_field = b[2]

    print first_field, second_field
}' \
> "all_haplotypes_taxid.txt"


# Keep only the taxa ID belonging to acceptable species.
# In this case, skip the first row (header) and extract the first column (TaxId).
tail --lines +2 "$acceptable_species_file" \
| cut --fields 1 \
| sort --numeric-sort \
| uniq \
> "acceptable_taxid.txt"


# Match and filter the haplotype file with the list of acceptable taxa IDs
awk \
'NR==FNR {nums[$1]; next} $2 in nums' \
"acceptable_taxid.txt" \
"all_haplotypes_taxid.txt" \
| awk \
--field-separator ' ' \
'{print $1, "taxid="$NF";"}' \
> "acceptable_haplotypes_taxid.txt"


# Extract only the list of acceptable haplotypes
cut \
--fields 1 \
--delimiter ' ' \
"acceptable_haplotypes_taxid.txt" \
> "acceptable_haplotypes.txt"


# Extract acceptable FASTA reads
zcat "$prefix.fasta.gz" \
| awk \
'/^>/{split($0, a, "."); print a[1] "." a[2]; next}1' \
| seqkit grep \
--pattern-file "acceptable_haplotypes.txt" \
--out-file "$prefix.temp.fasta"

awk \
--field-separator ' ' \
'NR==FNR {info[$1]=$2; next} /^>/ {id=substr($1, 2); print $1" "info[id]; next} 1' \
"acceptable_haplotypes.txt" \
"$prefix.temp.fasta" \
> "$prefix.filtered.fasta"

rm "$prefix.temp.fasta"

echo 'FASTA file filtered successfully'


# Download and prepare NCBI taxonomy dump
mkdir --parents "$taxonomy_dir"
cd "$taxonomy_dir"

wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'

tar \
--extract \
--overwrite \
--file \
'taxdump.tar.gz'

rm --force 'taxdump.tar.gz'

python3 \
"$blast_taxonomy_script" \
names.dmp \
nodes.dmp

echo "Taxonomy database ready"


# Make BLAST database
cd "$blast_out_dir"

# Map all haplotypes to their corresponding taxa ID, removing unwanted characters
sed 's/\(taxid=\|;\)//g' \
"$midori_dir/acceptable_haplotypes_taxid.txt" \
> "taxid_map.txt"

# Generate the BLAST database from the filtered FASTA file and taxonomic map
makeblastdb \
-in "$midori_dir/$prefix.filtered.fasta" \
-out "$blast_out_dir/$prefix" \
-dbtype nucl \
-parse_seqids \
-taxid_map "taxid_map.txt" \
-blastdb_version 5 \
-title "$prefix"

echo "BLAST databases generated successfully"


# Print some debug info (optional)
echo -e "\n\n\n---------- JOB INFO ----------\n\n"
qstat -f -1 "$PBS_JOBID"
