# Mitochondrial Metagenomics Tutorial
## Required Software
Unix OS or Unix Emulator<br/>Fastqc<br/>trimmomatic<br/>Spades<br/>BLASTn and BLASTx<br/>Quast<br/>MITOS<br/>python

## Install in your OS basic environment or create a custom environment for 
conda create --name ENVIRONMENT_NAME python=3.9 ipykernel -y

## Activate Genomics environment 
conda activate ENVIRONMENT_NAME

## Create  two column table of tab delimited lists identifying matching well and sample IDs
### Should look something like this:
A1	1EA<br/>B1	1EB<br/>C1	1EC<br/>D1	1ED

## Attach sampleIDs to their corresponding sequence file
### Fill the list file with tab separated pairs of well IDs to sample IDs
for file in *.fastq
<br/>do
<br/>wellID=$(awk -F"_" '{print $1}' <(echo $file))
<br/>sampleID=$(grep "$wellID" list | cut -f2) 
<br/>mv $file "$sampleID"_$file
<br/>done
<br/>cd "$working_directory"
<br/>done


## Quality Check
mkdir fastqc_raw-reads
<br/>fastqc *_R1_*.fastq *_R2_*.fastq -o fastqc_raw-reads

## Trim Reads and Quality Filter
mkdir trimmed-reads
<br/>trimmomatic PE -threads 8 SampleX_R1_001.fastq SampleX_R2_001.fastq \
    trimmed-reads/forward_reads_output\
   trimmed-reads/ reverse_reads_output\
    ILLUMINACLIP:adapters.fa:2:30:10:8:true\
    LEADING:3 TRAILING:3\
    SLIDINGWINDOW:4:10 MINLEN:36

## Assemble into contigs
metaspades.py -1 *R1_001.fastq -2 *R2_001.fastq -o metaspades_assembly

## MEGA BLAST contigs
blastn \
<br/>-task megablast \
<br/>-query assembled_contigs.fasta \
<br/>-db /home/genome/databases/nt \
<br/>-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
<br/>-culling_limit 5 \
<br/>-num_threads 48 \
<br/>-evalue 1e-10 \
<br/>-out SampleX.vs.nt.cul5.1e10.megablast.out

## Extract and ID ribosomal sequences
### Nuclear Ribosomal Sequences
barrnap --kingdom euk --threads 12 --outseq rRNA-seqs.fasta contigs.fasta
### Mitochondrial Ribosomal Sequences
barrnap --kingdom mito --threads 12 --outseq rRNA-seqs.fasta contigs.fasta

## Assess Remaining Contigs
quast.py contigs.fasta -o quast_results

## Identify the most promising Nematode mitochondrial sequences 
### BLASTn against the Nema-mtDB
blastn \
<br/>-query mito-contigs.fasta \
<br/>-subject Nema-mt-DB_v1.0.fasta \
<br/>-oufmt "6 qseqid sseqid pident qcovs length qstart qend sstart send gaps bitscore evalue" \
<br/>-out blastn_output.tsv

### Extract contig IDs that are most likely Nematode using BLASTn hits that match your chosen cutoffs
awk -v cutoff1="98" '$3 >= cutoff1 {print}' blastn_output.tsv | awk -v cutoff2="80" '$4 >= cutoff2 {print}'| cut -f2 -d":" | head -n1 | awk -F"|" '{print $(NF-1),$NF}'| cut -f1 >nematode-mito-contigids.list

### Extract the top taxonomic BLASTn hits 
while read -r line
<br/>do
<br/>grep -w "$line" blastn_output.tsv | head -n1
<br/>done<nematode-mito-contigids.list | cut -f2

### Extract hits that match an expected taxon
while read -r line
<br/>do
<br/>grep -w "$line" blastn_output.tsv
<br/>done<species.list | cut -f2

