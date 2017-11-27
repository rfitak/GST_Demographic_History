# SNP identification and PSMC reconstruction
In this section, SNPs are identified using a combination of [SAMTOOLS v1.3](https://github.com/samtools/samtools) and [BCFTOOLS v1.3](http://samtools.github.io/bcftools/).  After several filtering steps, the PSMC analysis is performed and results plotted.

## Step 1:  Call SNPs
```bash
# Call and Filter SNPs
samtools1.3 mpileup \
   -C 50 \
   -q 20 \
   -Q 20 \
   -d 200 \
   -u \
   -f C_mydas.fa \
   C_mydas.merged.sorted.bam | \
   bcftools1.3 call \
      -m -v | \
   bcftools1.3 filter \
      -g10 -G10 -O v | \
   bcftools1.3 filter \
      -sFAIL -e'%QUAL<20 || \
      INFO/DP<=15 || \
      INFO/DP>=88 || \
      INFO/MQ<=30 || \
      INFO/MQB<1e-20 || \
      INFO/RPB<0.0001 || \
      INFO/BQB<0.0001 || \
      (INFO/DP4[1]+INFO/DP4[2]<=2) || \
      (INFO/DP4[3]+INFO/DP4[4]<=2)' \
      -O v > C_mydas.filtered.vcf      

# Remove INDELS, Mitochondrial SNPs, and retain only "PASS" SNPs
grep "^#" C_mydas.filtered.vcf | \
   cat - <(grep -v "^#" ../C_mydas.filtered.vcf | grep "PASS" | grep -v "INDEL" | grep -v "NC_000886.1") \
   > input.vcf
   # Dont forget to make a new reference fasta file and index
   # using the genome without the mitochondrial contig
```
Description of parameters:
- samtools1.3 mpileup  ::  create pileup file
    - -C 50  :: adjust mapping quality; recommended:50
    - -q 20  :: skip alignments with mapQ smaller than 20
    - -Q 20  :: skip bases with baseQ/BAQ smaller than 20
    - -d 200  :: max per-BAM depth (200x coverage); avoids excessive memory usage
    - -u  :: uncompressed BAM output
    - -f C_mydas.fa  :: indexed reference genome file
- bcftools1.3 call  :: SNP/indel variant calling from VCF/BCF. To be used in conjunction with samtools mpileup
    - -m  :: alternative model for multiallelic and rare-variant calling (conflicts with -c)
    - -v  :: output variant sites only
- bcftools1.3 filter  :: Apply fixed-threshold filters to VCF file
    - -g10  :: filter SNPs within 10 base pairs of an indel
    - -G10  :: filter clusters of indels separated by 10 or fewer base pairs allowing only one to pass
    - -O v  :: output type v: uncompressed VCF
    - -sFAIL :: annotate FILTER column with "FAIL"
    - -e'%QUAL<20  :: exclude sites for which the expression is true (see man page for details), SNP Qaulity < 20
    - INFO/DP<=15  :: read depth <= 15X (1/3 genome-wide average)
      INFO/DP>=88  :: read depth >= 88X (2X genome-wide average)
      INFO/MQ<=30  :: mean mapping quality of SNP reads <=30
      INFO/MQB<1e-20  :: mapping quality strand bias p value less than 10^-20
      INFO/RPB<0.0001  :: read position bias p value < 0.0001
      INFO/BQB<0.0001  :: base quality bias p value < 0.001
      (INFO/DP4[1]+INFO/DP4[2]<=2)  :: at least 2 reads with the reference base
      (INFO/DP4[3]+INFO/DP4[4]<=2)'  :: at least 2 reads with the alternate base
      
      
      
## Step 2:  Create consensus FASTA file
Here we apply the ambiguous base codes to all SNPs in the reference genome and and convert the output into a fasta file used by PSMC.

```bash
# Make Consensus fasta file
bcftools1.3 consensus \
   -f C_mydas.fa \
   -i input.vcf | \
   gzip > C_mydas.filtered.fa.gz
```
Description of parameters:
- bcftools1.3 consensus  :: Create consensus sequence by applying VCF variants to a reference fasta file
    -f C_mydas.fa  :: indexed reference genoe file
    -i  :: output variants in the form of IUPAC ambiguity codes



## Step 3:  Run PSMC
This step will generate the PSMC-specific input file format, run PSMC, and also run PSMC on 100 bootstrap replicates.

```bash

```



