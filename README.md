# Aqua-eDNA-metabarcoding-pipeline

## I. Preliminary steps

### 1. Taxonomy

- **Website**: [SeaLifeBase](https://www.sealifebase.org/search.php)

Having access to a curated list of species we might encounter in specific areas can tremendously increase the accuracy of taxonomic assignment when targeting short mitochondrial regions.
For example, the Atlantic cod and Pacific cod share more than 99% of their mitochondrial DNA, and it might be impossible to distinguish them from the 12S subunit alone.
By creating a list of acceptable taxa, we can filter out the noise caused by closely genetically related species residing far away from our collection area.

We should aim to create a table of acceptable species containing at least two variables: the NCBI taxonomy ID and the corresponding scientific name.
The common name can also be helpful for data visualisation, together with other project-specific features. This is a small example:

| **TaxId** | **Scientific_name**   | **Common_name**        | **Category**     |
|-----------|-----------------------|------------------------|------------------|
| `7748`	   | Lampetra fluviatilis	 | European river lamprey | 	Fishes          |
| `7757`	   | Petromyzon marinus    | 	Sea lamprey           | 	Fishes          |
| `7769`	   | Myxine glutinosa      | 	Atlantic hagfish      | 	Fishes          |
| `28701`   | Fratercula arctica    | Atlantic puffin        | Birds            |
| `30455`   | Fulmarus glacialis    | Northern fulmar        | Birds            |
| `9707`    | Odobenus rosmarus     | Walrus                 | Marine mammals   |
| `9711`    | Halichoerus grypus    | Gray seal              | Marine mammals   |
| `9031`    | Gallus gallus         | Chicken                | Contamination    |
| `9606`    | Homo sapiens          | Human                  | Contamination    |
| `30501`   | Carcharias taurus     | Sand tiger shark       | Positive control |

As a starting point, we used [SeaLifeBase](https://www.sealifebase.org/search.php) to get a list of common marine species residing in Denmark.
We also included additional elements relevant to our study.
Make sure to add natural contaminants, with genetic material stemming from human activities and sewage wastewater, as an example.
We encountered a high amount of reads assigned to Humans, Chickens, Pigs, Domestic cattle, and Sheep.

For reference, you can find the full table we used for our study [by clicking here](sample_data/acceptable_species.tsv).

#### Core commands

```bash
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
```

---------------------------------------------------

## II. Metabarcoding and OTU formation

### 1. Merge paired-end reads

- **Tool**: [PEAR](https://cme.h-its.org/exelixis/web/software/pear/doc.html)
- **Full script**: [1.pear.pbs.sh](scripts/1.pear.pbs.sh)

PEAR is used to merge paired-end reads contained in separate FASTQ files. Unassembled reads are discarded.

#### Core commands

```bash
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
```

---------------------------------------------------

### 2. Check quality of assembled reads

- **Tool**: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- **Full script**: [2.fastqc.pbs.sh](scripts/2.fastqc.pbs.sh)

The assembled FASTQ file is passed to FastQC to perform semi-quantitative quality control of the merged reads.

#### Core commands

```bash
fastqc \
--outdir "$out_dir" \
--threads 10 \
--quiet \
--dir "$out_dir" \
"$fastq_dir/LibName.assembled.fastq"
```

---------------------------------------------------

### 3. Metabarcoding pipeline

- **Tools**: [Cutadapt](https://cutadapt.readthedocs.io/en/stable/), [OBITools](https://pythonhosted.org/OBITools/index.html), [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/)
- **Full script**: [3.metabarcoding.pbs.sh](scripts/3.metabarcoding.pbs.sh)
- **Additional files required**: adapter files, BLAST database, taxonomy information, [`taxonomy_assignment_BLAST.py`](scripts/taxonomy_assignment_BLAST.py)

This is the core pipeline for metabarcoding and taxonomic assignment. It can be split into the following steps:

1. **`cutadapt`**: Individual reads are split by sample tags. The `sample` and `tag_matches` parameters are added to all headers of the FASTA file, reporting the sample name and matched sequences.
2. **`cutadapt`**: Individual reads are split by primers. The `primer` and `primer_matches` parameters are added to all headers of the FASTA file, reporting the primer name and matched sequences.
3. **`obiuniq`**: Identical reads are merged to compress the FASTA file and speed up the computation.
4. **`obistat`**: Compute basic statistics for read counts (optional).
5. **`obigrep`**: Remove read groups with less than 2 matches or shorter than 20 nucleotides.
6. **`obiclean`**: Remove spurious reads and potential errors caused by PCR or sequencing.
7. **`obiannotate`**: Prune unused metadata from read headers (optional).
8. **`obitab`**: Extract metadata from FASTA headers and export it to a TSV file for downstream analysis.
9. **`blastn`**: Perform a first round of taxonomic assignments for all OTUs.
10. **`taxonomy_assignment_BLAST.py`**: Refine BLAST results by applying a consensus algorithm to the raw data from `blastn`.

#### Core commands

<details><summary>Show code</summary>

```bash
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

# [...]

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

# [...]

obiuniq \
--merge sample \
--without-progress-bar \
"2_demultiplex/$library/$library.trimmed.fastq" \
> "$out_dir/$library.unique.fasta"

# [...]

obistat \
--category-attribute count \
--without-progress-bar \
"$out_dir/$library.unique.fasta" \
| sort --numeric-sort --key 1 \
--output "$out_dir/$library.unique.stats.txt"

# [...]

obigrep \
--without-progress-bar \
--lmin 20 \
--predicat 'count >= 2' \
"3_obiuniq/$library.unique.fasta" \
> "$out_dir/$library.pre-filtered.fasta"

# [...]

obiclean \
--without-progress-bar \
--sample merged_sample \
--ratio 0.05 \
--head \
"$out_dir/$library.pre-filtered.fasta" \
> "$out_dir/$library.clean.fasta"

# [...]

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

# [...]

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

# [...]

python3 \
"$blast_taxonomy_script" \
--hits_to_consider 5 \
--length_percentage 0.97 \
--percent_sway 0.5 \
--ncbi_nt \
--blast_database IGNORE \
--blast_file "$out_dir/$library.raw.blast.tsv" \
--output_dir "$out_dir" \
"$fasta_file" \
"$taxonomy_file"
```

</details>

## III. Data analysis and plots

