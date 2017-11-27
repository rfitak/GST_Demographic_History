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

### Step 1: Trim sequences
Next, we trim low-quality bases and reads using [TRIMMOMATIC v0.35](http://www.usadellab.org/cms/?page=trimmomatic).

```bash
# Trim reads
trimmomatic \
   PE \
   -threads 8 \
   input_F.fastq.gz \
   input_R.fastq.gz \
   output_F.fastq.gz \
   output_F.SE.fastq.gz \
   output_R.fastq.gz \
   output_R.SE.fastq.gz \
   ILLUMINACLIP:all.fa:2:30:7 \
   LEADING:20 \
   TRAILING:20 \
   SLIDINGWINDOW:4:20 \
   AVGQUAL:30 \
   MINLEN:50

```
Desciption of the parameters:
- PE :: the input data are paired-end reads
- -threads 8 :: use four CPUs (threads)
- input_F.fastq.gz :: the name of the forward reads file
- input_R.fastq.gz :: the name of the reverse reads file
- name of the trimmed and paired forward reads output file
- name of the trimmed and unpaired forward reads output file
- name of the trimmed and paired reverse reads output file
- name of the trimmed and unpaired reverse reads output file
- ILLUMINACLIP:all.fa:2:30:7 :: fasta file of adapter sequences to trim, see file [all.fa](./Data/all.fa)
  * 2 :: seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
  * 30 :: palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
  * 7 :: simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.
 - LEADING:20 :: trim bases from start of the read with a Q < 20
 - TRAILING:20 :: trim bases from end of the read with a Q < 20
 - SLIDINGWINDOW:4:20 :: using a 4-base window, remove the last base if average Q < 20
 - AVGQUAL:30 :: remove the entire read if the average quality is < 30
 - MINLEN:50 :: remove reads less than 50 bases after quality trimming
