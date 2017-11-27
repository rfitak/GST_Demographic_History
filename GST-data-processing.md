# Downloading and Processing Raw Sequence Data


### Step 1: Download sequence data
Here we download only the raw, short-insert, paired-end sequencing data. Short insert sequences are recommended for SNP identification and thus any downstream SNP-related analyses.  Long-insert sequences are primarily for scaffolding and structural analyses.  To download the data from the Sequence Read Archive (SRA database), we use the NCBI [SRATOOLKIT](https://github.com/ncbi/sra-tools).

The Run Table describing all the SRA metadata for the study can be downloaded here [SAMN02981410 Run Table](http://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN02981410).  It also available here as [SraRunTable.txt](./Data/SraRunTable.txt).

```bash
# Get list of Run IDs
sed '1d' SraRunTable.txt | cut -f7 | head -8 > SRR.list
```

The file, [SRR.list](./Data/SRR.list], contains the run IDs for the 8 sequencing libraries to be downloaded. Now they can be downloaded using the commands below.

```bash
while read line
   do
   fastq-dump --split-files -Z -gzip $line > $line.fq.gz
   echo "Finished $line"
done < SRR.list
```
Description of parameters:
    - --split-files :: separate forward and reverse read pairs into separate files
    - -Z            :: write to standard out
    - -gzip         :: compress output using gzip



```bash
1.	Trim reads using trimmomatic v0.35:
a.	trimmomatic PE -threads 8 input_F.fastq.gz input_R.fastq.gz output_F.fastq.gz output_F.SE.fastq.gz output_R.fastq.gz output_R.SE.fastq.gz ILLUMINACLIP:all.fa:2:30:7 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 AVGQUAL:30 MINLEN:50

```
