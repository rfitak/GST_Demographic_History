# SNP identification and PSMC reconstruction

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
