Scripts:

10_filter_table_by_threshold.py
Based on the results of the MASH algorithm, filter out similar genomes.

11_run_prokka_annotation.sh

This script runs Prokka in batch mode to annotate genome files (.fna). It generates amino acid sequences (.faa), nucleotide sequences (.ffn), annotation tables (.tsv), and reformatted genome files (.fna with line length changed from 80 to 60). Parallel processing is used to enhance efficiency.

12_fetch_cds_dna_seq_bacteria.R

This R script extracts CDS (coding sequences) DNA sequences from Prokka annotation files (.tsv) and saves them as FASTA files. The script checks if the output file already exists to avoid redundant processing and logs any errors encountered during execution.


**Filter genes with a size greater than 10k

13_extract_long_genes.sh
This script extracts gene sequences longer than 10,000 bases. It uses "seqkit" to convert FASTA files to tabular format, extracts the length of each gene, and saves the IDs of genes longer than 10,000 bases to a text file. The script supports parallel processing to enhance efficiency.

14_fetch_specifi_gene_fq.py
This Python script extracts specific gene sequences from FASTA files. It reads the input FASTA file and a list of gene IDs, outputting sequences that match the criteria to one file and those that do not to another.

15_fetch_specifi_gene_fq.sh
This script extracts specific gene sequences from FASTA files. It uses the Python script "fetch_specifi_gene_fq.py" to process each gene file, categorizing the sequences into two groups based on a provided gene ID list: sequences that match the criteria and those that do not. The script supports parallel processing to enhance efficiency.

**Remove redundancies using Cd-hit

16_merge_fa.sh
This script merges multiple FASTA files into several larger files. It processes all FASTA files in the input directory in batches, with each batch containing 1000 files. The merged files are saved to the specified output directory. The script uses a temporary file to store the list of files for each batch and deletes the temporary file after processing.

17_cdhit_cluster_analysis.sh
This script is used to batch-run the cd-hit-est tool for cluster analysis of sequence files in FASTA format.

18_merge2fasta.py
This script splits merged FASTA sequence files into genome-specific FASTA files. It loads merged FASTA files from a specified directory, reads gene name text files to obtain genome-specific gene IDs, matches sequences based on these IDs, and generates dedicated FASTA files for each genome. It supports batch processing with command-line parameters for input/output paths.

19_rebuild_fasta.sh
This script executes the merge2fasta.py Python script. It activates the specified conda environment, defines paths for input (merged FASTA files, gene name files), output directory, and the Python script itself, then runs the Python script with necessary parameters to perform sequence splitting.

**Split CDS
20_split_cds_to_pseudo_reads.py
Each CDS sequence is sliced into "pseudo reads" with a 150-bp non-overlapping sliding window, and the FASTQ (including simulated quality) and the statistics table of the number of fragments per gene are output.
21_run_split_cds_to_pseudo_reads.sh
Start the "split_cds_to_pseudo_reads.py" script and import the relevant parameters.

**Bowtie2
```bash
bowtie2-build -f --large-index --bmax 6635772616 --dcv 4096 --threads 28 bacteria.fna bacteria
```

22_map_cds_pseudoreads_bt2.sh
Map CDS pseudo-reads (FASTQ) to the bacterial Bowtie2 index in end-to-end mode; keep uniquely aligned reads (AS:i present, XS:i absent) and generate SAM, aligned/unaligned FASTQ plus logs.

23_fetch_specific_gene_name.py
24_start_fetch_specific_gene_name.sh

In the comparison results, only the CDS genes with "all 150-bp fragments are uniquely aligned" are retained, and the list of their gene IDs is output.

**Uniref Annotator

```bash
#Generate index file by diamond
diamond makedb --in /s3/SHARE/woodman/Prokka2/data/uniref90.fasta --db /s3/SHARE/woodman/Prokka2/dmnd_db/uniref90.dmnd
diamond makedb --in /s3/SHARE/woodman/Prokka2/data/uniref50.fasta --db /s3/SHARE/woodman/Prokka2/dmnd_db/uniref50.dmnd
```

25_batch_diamond_uniref.sh
26_start_mutiple.sh
Takes a sample list, launches 2 parallel jobs (7 CPUs each) to DIAMOND-search 40 k FAA files against UniRef90 first and UniRef50 next. Outputs are saved under result/4w/<sample>/; finished samples are skipped automatically.

27_extract_genes_to_individual_fasta.py
28_run_extract_genes_to_individual_fasta.sh
According to the list of gene IDs, split the sequences in the input FASTA file into individual {gene_id}.fasta files and generate a completion log.

**Primer Design
```bash
conda install -c bioconda primer3-py=0.6.1
conda install bioconda::pysam
conda install cctbx202105::biopytho
conda install pandas=1.5.3 numpy=1.21.2
```
29_primer_design.py
30_run_primer_design.sh

For each input gene sequence, call Primer3 to batch design qPCR primers (with product sizes ranging from 100 to 500 bp), and output the CSV results and log files, skipping sequences shorter than 100 bp.

31_consolidate_primers.py

Scan the virus primer result directory, horizontally merge all single-gene CSV files with a size greater than 3 bytes, add the GENOME_ID / GENE_ID columns, and output the final table named "primer_bank_virus.csv"

32_score_primers.py
33_run_score_primers.sh

Read the primer CSV file, calculate ΔG, GC%, Tm, hairpin/self-complementary penalty for each pair, and output the comprehensive score and intermediate indicators.


34_primer_merge.py
35_primer_mergeV2.py
Recursively read all primer CSV files in the "ann_uniref" directory, batch append and merge them into a single "primer_bank_bacteria_4w_v1.csv" file, and automatically complete the unique indexes for KINGDOM and PRIMER_PAIR_X.


36_gtf2json.py
37_process_gff.py
Batch convert the "gene" rows in GTF into JSONL format: Convert each row to a dictionary → Add a file prefix → Write to individual files first and then aggregate into gtf.jsonl. Taking fungi as an example.

38_assembly.json.merge.py
Convert the two TSV files (NCBI bacteria assembly summary (historical + current)) into GTSS (Genome-TAXID-Species-Strain) JSONL format. First, write them into separate files, then merge them into a non-redundant master database.

39-43_primer_merge.py 
Update form fields and data(V2 to v7)

**nt
44_map_nt_bowtie2_unique.sh
End-to-end map 40 k CDS pseudo-sequences against the NCBI nt partition index with bowtie2; keep only uniquely aligned reads (AS:i present, XS:i absent) and produce SAM, aligned/unaligned FASTA plus logs, accelerated by 50 parallel jobs.
```bash
find nt_A_map -type f -name "*.sam" -exec cp {} ./nt_A/  \;

find nt_B_map -type f -name "*.sam" -exec cp {} ./nt_B/  \;
```

45_merge_sam.py
Merge the comparison results of the nt A and B partitions by genomic prefix, remove duplicates, and eliminate the rows containing "XS:i", and output the unique aligned SAM file.

46_ncbi_refseq_batch_dl_md5.py
Conduct multi-threaded batch download of NCBI RefSeq genome sequence reports, extract and verify MD5 in real time, and record the failed IDs and verification errors separately.

47_sam_to_region_primer_index.py
48_sam_to_region_primer_index_v2.py
Scan the "merged_sam" directory, extract the region and primer name for each alignment record, output a JSON index for each genome and summarize them into a TSV/JSON/PKL overall table.

49_merge_seq_reports_index.py
Iterate through all *_seq_report.jsonl files in the sequence_reports directory, perform deduplication based on genbankAccession, and summarize them into the assembly/refseq/moleculeType index. Output the JSON/PKL summary table and the single genome JSON cache.

**Generate Download File
50_parallel_tar_taxid_dirs.sh
Parallelly package all subdirectories under "temp_fna_taxid" into "${subdir}.fna.tar.gz"

51_generate_6format_taxid_packages.py
Split and package 6 types of download files according to TAXID: csv (already split), faa/fa gene groups and CDS sequences, gff→bed. Real-time alerts for missing records.

```bash
source /home/caisongbo/anaconda3/bin/activate primer

mkdir -p gb

cd genome_browser

for f in *.fna; do
    samtools faidx "$f"
    prefix=${f%.fna}

    if [ -d "../gb/${prefix}" ]; then
        echo "Directory ../gb/${prefix} already exists, skipping."
        continue
    fi

    mkdir -p ../gb/"${prefix}"
    mv "${prefix}.fna" ../gb/"${prefix}"
    mv "${prefix}.fna.fai" ../gb/"${prefix}" 2>/dev/null
    mv "${prefix}.gff" ../gb/"${prefix}" 2>/dev/null

done

mate_virus/
├── bed_taxid
├── csv_taxid
├── faa_taxid
├── fa_taxid
├── fna_taxid
└── gff_taxid

csv_taxid/
1000373.csv.gz
1000646.csv.gz
1001080.csv.gz

gff_taxid/
1000373.gff.tar.gz
1000646.gff.tar.gz
1001080.gff.tar.gz
$ tar -tzf gff_taxid/1001080.gff.tar.gz
GCF_002830685.1.gff

fna_taxid/ 
1000373.fna.tar.gz
1000646.fna.tar.gz
1001080.fna.tar.gz
$ tar -tzf fna_taxid/1002921.fna.tar.gz
GCF_000863425.1.fna

bed_taxid/
1000373.bed.tar.gz
1000646.bed.tar.gz
1001080.bed.tar.gz
bed_taxid
GCF_001551505.1.bed

fa_taxid/ 
1000373.fa.tar.gz
1000646.fa.tar.gz
1001080.fa.tar.gz
$ tar -tzf fa_taxid/1003177.fa.tar.gz
GCF_001551505.1.fasta

faa_taxid/
1000373.faa.tar.gz
1000646.faa.tar.gz
1001080.faa.tar.gz
$ tar -tzf faa_taxid/1000646.faa.tar.gz
GCF_013086425.1.faa
```
