# Downloading and Processing Raw Sequence Data


## Step 1: Download sequence data
Here we download only the raw, short-insert, paired-end sequencing data. Short insert sequences are recommended for SNP identification and thus any downstream SNP-related analyses.  Long-insert sequences are primarily for scaffolding and structural analyses.  To download the data from the Sequence Read Archive (SRA database), we use the NCBI [SRATOOLKIT v2.8.2](https://github.com/ncbi/sra-tools).

The Run Table describing all the SRA metadata for the study can be downloaded here [SAMN02981410 Run Table](http://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN02981410).  It also available here as [SraRunTable.txt](./Data/SraRunTable.txt).

```bash
# Get list of Run IDs
sed '1d' SraRunTable.txt | cut -f7 | head -8 > SRR.list
```

The file, [SRR.list](./Data/SRR.list], contains the run IDs for the 8 sequencing libraries to be downloaded. Now they can be downloaded using the commands below.

```bash
while read SRA
   do
   fastq-dump --split-files -gzip $SRA
   echo "Finished $SRA"
done < SRR.list
```
Description of parameters:
- --split-files :: separate forward and reverse read pairs into separate files
- -gzip         :: compress output using gzip



## Step 2: Trim sequences
Next, we trim low-quality bases and reads using [TRIMMOMATIC v0.35](http://www.usadellab.org/cms/?page=trimmomatic).  The code loops through each of the 8 pairs of SRA files (16 in total) separately.

```bash
# Trim reads, loop through each SRA file
while read SRA
   do
   trimmomatic \
      PE \
      -threads 8 \
      ${SRA}_1.fastq.gz \
      ${SRA}_2.fastq.gz \
      ${SRA}_F.trimmed.fastq.gz \
      ${SRA}_F.SE.trimmed.fastq.gz \
      ${SRA}_R.trimmed.fastq.gz \
      ${SRA}_R.SE.trimmed.fastq.gz \
      ILLUMINACLIP:all.fa:2:30:7 \
      LEADING:20 \
      TRAILING:20 \
      SLIDINGWINDOW:4:20 \
      AVGQUAL:30 \
      MINLEN:50
   done < SRR.list
```
Desciption of the parameters:
- PE :: the input data are paired-end reads
- -threads 8 :: use four CPUs (threads)
- ${SRA}\_1.fastq.gz :: the name of the forward reads file
- ${SRA}\_2.fastq.gz :: the name of the reverse reads file
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



## Step 3: Download and index reference genome
The *C. mydas* genome is downloaded from GenBank and then indexed using [BWA v0.7.12](http://bio-bwa.sourceforge.net).

```bash
# Download reference genome
wget https://ftp.ncbi.nih.gov/genomes/Chelonia_mydas/CHR_Un/8469_ref_CheMyd_1.0_chrUn.fa.gz
gunzip 8469_ref_CheMyd_1.0_chrUn.fa.gz
mv 8469_ref_CheMyd_1.0_chrUn.fa C_mydas.fa

# Build Index
bwa index \
   -a bwtsw \
   C_mydas.fa
```



### Step 4: Map trimmed reads to reference genome

```bash
# Define a read group in header
rg="@RG\tID:Cmydas\tSM:Cmydas"

# Begin mapping, looping through each SRA file pair of trimmed reads
while read SRA
   do
   bwa mem \
      -t 8 \
      -R $rg \
      -M \
      C_mydas.fa \
      ${SRA}_F.trimmed.fastq.gz \
      ${SRA}_R.trimmed.fastq.gz | \
      samtools1.3 view \
         -b \
         -q 20 \
         -f 0x0002 \
         -F 0x0004 \
         -F 0x0008 \
         -T C_mydas.fa \
         - | \
      samtools1.3 sort -T ${SRA}.tmp \
         - | \
      samtools rmdup \
         - ${SRA}.sorted.bam
   
   # Index the BAM file
   samtools1.3 index ${SRA}.sorted.bam

   # Get mapping stats from the BAM file
   samtools1.3 stats ${SRA}.sorted.bam > ${SRA}.bamstats
   
done < SRR.list
```
Desciption of the parameters:
- -t 8  :: 8 threads
- -R $rg  :: add the read group line to the BAM header
- -M  :: mark shorter split hits as secondary
- C_mydas.fa  :: indexed reference genome file
- samtools1.3 view  :: use samtools v1.3 to filter aligned reads
    - -b  :: xxxx
    - -q 20  :: xxxx
    - -f 0x0002  :: xxxx
    - -F 0x0004  :: xxxx
    - -F 0x0008  :: xxxx
    - -T C_mydas.fa  :: indexed reference genome file
- samtools1.3 sort  :: use samtools v1.3 to sort the BAM file
    - -T ${SRA}.tmp  :: store temporary sorted files here
samtools rmdup  :: use the older samtools v0.1.19 to remove duplicate reads

