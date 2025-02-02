# eDNA metabarcoding pipeline

This repository collects scripts and instructions useful for performing metabarcoding analyses from environmental DNA (eDNA).
Paired-end FASTQ reads will be clustered in OTUs and assigned to the most likely taxonomic group.
The data can then be plotted and studied to infer more specific trends.

For each step, you will find an overview, a link to a full Bash or R script, and, if relevant, examples of what you will obtain.

### Use case

This pipeline was developed in collaboration with the [National Institute of Aquatic Resources](https://www.aqua.dtu.dk/english/) (DTU Aqua).
It was used as part of the Master's thesis project: "[Monitoring offshore biodiversity through the collection of eDNA using an autonomous ecogenomic sensor](https://rawcdn.githack.com/ManuelArcieri/Aqua-eDNA-metabarcoding-pipeline/main/sample_data/Thesis.pdf)", by Manuel Arcieri.

The sample data contains information from a campaign conducted by DTU Aqua in the North Sea (2019).
Water samples were collected from several sites and analysed for taxonomic assignment of bony fish. The [thesis](https://rawcdn.githack.com/ManuelArcieri/Aqua-eDNA-metabarcoding-pipeline/main/sample_data/Thesis.pdf) provides more details.

## I. Preliminary steps

Before you can analyse your samples, you have to define a reference database for taxonomic assignment.
Custom databases can be more time-consuming to create, but they can provide more pertinent results.

### 1. Taxonomy

- **Website**: [SeaLifeBase](https://www.sealifebase.org/search.php)

Having access to a curated list of species you might encounter in specific areas can tremendously increase the accuracy of taxonomic assignment when targeting short mitochondrial regions.
For example, the Atlantic cod and Pacific cod share more than 99% of their mitochondrial DNA, and it might be impossible to distinguish them from the 12S subunit alone.
By creating a list of acceptable taxa, you can filter out the noise caused by closely genetically related species residing far away from our collection area.

You should aim to create a table of acceptable species containing at least two variables: the NCBI taxonomy ID and the corresponding scientific name.
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

#### Iterative refinement

As part of the metabarcoding pipeline, you will run BLAST on your sequences. You can inspect the raw results to select high-frequency OTUs that couldn't be assigned to any taxonomic group.
Here's an example:

| **OTU** | **Count** | **Prediction_rank** | **Prediction_s_name**  | **Prediction_c_name** | **...** |
|---------|-----------|---------------------|------------------------|-----------------------|---------|
| A       | 10260	    | Species	            | Clupea harengus        | 	Atlantic herring     | ...     |
| B       | 9461	     | Species	            | Clupea harengus        | 	Atlantic herring     | ...     |
| C       | 8342	     | Species	            | Bos taurus             | 	Domestic cattle      | ...     |
| **D**   | **6104**	 | **NA**	             | **NA**                 | 	**NA**               | ...     |
| E       | 5719	     | Species	            | Bos taurus             | 	Domestic cattle      | ...     |
| F       | 5511	     | Species	            | Pleuronectes platessa	 | European plaice       | ...     |

As you can see, OTU `D` wasn't assigned to any taxa.
Usually, this is due to a missing species or haplotype in the list of acceptable elements or reference database.
You can manually extract the DNA sequence from the FASTA file matching the unclassified OTU.
Then, you can search it online using the [full NCBI database](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch) (`blastn`) to determine the most likely species it belongs to.
You can add it to your list and re-run the taxonomic assignment of all your OTUs.
This procedure can be repeated until you end up with satisfying results.

---------------------------------------------------

### 2. Reference database

- **Website**: [MIDORI2](https://www.reference-midori.info/index.html)

A reference database contains known DNA sequences from various organisms.
It is used for taxonomic assignment by comparing an unknown sequence with all those found in the database, producing similarity scores for the top matches.
Many reference databases may have multiple haplotypes associated with a single species, potentially increasing the accuracy of predictions.

The choice of a reference database depends on the targeted DNA regions and the scope of the project.
In our case, we opted for [MIDORI2](https://www.reference-midori.info/index.html).
It subsets the NCBI GenBank database, providing DNA and amino acid sequences for eukaryotic mitochondrial DNA regions.
We targeted the 16S region for marine mammals, and the 12S region for bony fish and elasmobranchs.

---------------------------------------------------

### 3. Creating a BLAST reference database (using MIDORI2)

- **Tools**: [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.makeblastdb_application_opt/), [SeqKit](https://bioinf.shenwei.me/seqkit/)
- **Full script**: [make-blast-database.pbs.sh](scripts/make-blast-database.pbs.sh)
- **Additional files required**: FASTA reference sequences, list of acceptable species, [`genbank_nodes_and_names_to_taxonomy.py`](scripts/genbank_nodes_and_names_to_taxonomy.py)

Once you have a list of acceptable species (e.g., [`acceptable_species.tsv`](sample_data/acceptable_species.tsv)) and a FASTA reference database (e.g., [MIDORI2](https://www.reference-midori.info/index.html) in QIIME format), you are ready to create
a custom BLAST database that can be used for taxonomic assignment.

All the steps are merged into a single script, [make-blast-database.pbs.sh](scripts/make-blast-database.pbs.sh), and are summarised below:

1. Extract all haplotypes and taxa IDs contained in the reference database.
2. Extract all acceptable taxa IDs from our TSV file.
3. Filter out undesired species/haplotypes.
4. Create a filtered FASTA file containing all acceptable haplotypes.
5. Download the latest NCBI taxonomy dump and prepare it for later usage.
6. Make a BLAST database from the filtered FASTA file and taxa IDs.

### 4. Preparing the adapter files for Cutadapt

- **Tool**: [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html#demultiplexing)

As part of the metabarcoding pipeline, you will run Cutadapt twice to trim the primers and demultiplex samples based on their tags.
A library may contain an arbitrary number of primers and tags. Thus, we need two files to list all of them.
This is a simple process. You will need two FASTA files where each primer/tag is listed as a single record.

The sequence can contain Cutadapt's special instructions.
For example, the form `AAA...XXX` can be used to specify a linked adapter, where your target sequence (`...`) is surrounded by a forward (`AAA`) and reverse (`XXX`) adapter.

In the FASTA files of the adapters, the header name must be unique and can be used to demultiplex different samples.

Primers typically do not require specific identifiers.
You can append a sequence number to differentiate them (e.g., [`LibDemo.primers.fasta`](sample_data/LibDemo.primers.fasta)):

```text
>primer1
GTCGGTAAAACTCGTGCCAGC...CAAACTGGGATTAGATACCCCACTATG
>primer2
CATAGTGGGGTATCTAATCCCAGTTTG...GCTGGCACGAGTTTTACCGAC
```

On the other hand, tags are typically used to discriminate samples or experimental conditions.
You can use more specific headers to simplify downstream data analysis (e.g., [`LibDemo.tags.fasta`](sample_data/LibDemo.tags.fasta)):

```text
>DanF_surface_1
ACCGAGA...TCTCGGT
>1600m_East_surface_1
CGAAATG...CATTTCG
>Ragnar_surface_1
TGAACGG...CCGTTCA
>Novana_shipwreck_surface_1
CAACCTT...AAGGTTG
```

Cutadapt will append the tag and primer IDs to the sequence headers during demultiplexing.
These fields are kept for subsequent analyses and can be inspected at any time. Here's an example:

```text
@EXAMPLE1 sample=1600m_East_surface_1; tag_matches=CGAAATG,CATTTCG; primer=primer1; primer_matches=GTCGGTAAAACTCGTGCCAGC,CAAACTGGGATTAGATAC
CACCACACGATTAACCCGAGCTAATGGAACT
+
F/FF6FAFF//F/FFFFF6FFFFFFAFF/F/
```

## II. Metabarcoding and OTU formation

Once you have defined your taxonomic target group, you can proceed with the metabarcoding analysis.

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
- **Additional files required**: adapter files, BLAST database, taxonomy dump, [`taxonomy_assignment_BLAST.py`](scripts/taxonomy_assignment_BLAST.py)

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

The computationally intensive part of the pipeline has now ended. The results can be further refined to improve readability and to plot the main features of the dataset.
The following steps can be run on a personal computer with [R](https://www.r-project.org) installed.

**(To be added...)**

## Acknowledgements

The bulk of the code was written by Manuel Arcieri, former Master's student at the Technical University of Denmark (DTU), as part of his thesis project.

I would also like to thank Sara Maggini (DTU Aqua) and Magnus Wulff Jacobsen (DTU Aqua) for their precious help.

If you have any questions, you can write to Manuel Arcieri at `manuel.arcieri@outlook.com`.

## Code usage

You are free to copy, modify, and distribute all the code from this repository as described in the [license](LICENSE).
