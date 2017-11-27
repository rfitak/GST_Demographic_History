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
# Generate psmcfa input file
fq2psmcfa C_mydas.filtered.fa.gz > C_mydas.filtered.psmcfa

# Split for bootstrapping
splitfa C_mydas.filtered.psmcfa > C_mydas.split.filtered.psmcfa

# Run PSMC on the genome
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o C_mydas.filtered.psmc C_mydas.filtered.psmcfa

# Run PSMC on the bootstrapped files
for b in {1..100}
   do
   psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o C_mydas.filtered.${b}.bs.psmc C_mydas.split.filtered.psmcfa
done

# Merge together full and bootstrapped files
cat \
   C_mydas.filtered.psmc \
   C_mydas.filtered.{1..100}.bs.psmc > \
   C_mydas.filtered.combined.psmc
```



## Step 4: Plot Results

```bash
psmc_plot.pl \
   -X50000000 \
   -p \
   -g42.8 \
   -R \
   -x1000 \
   -u1.2e-08 \
   test \
   C_mydas.filtered.combined.psmc
```


## Step 5:  Additional Plotting functions in R
This is extra code to further make specific plots frm the PSMC output, in addition to adding a sliding window of sea surface temperature data taken from [van de Wal et al. 2011 *Climate of the Past*](https://www.clim-past.net/7/1459/2011/).
```R
# Load libraries needed
library(scales)
library(zoo)

# Open PDF file to write to
pdf("psmc.plot.pdf", height = 7, width = 10)

# Adjust plot margins
par(mar = c(4.1, 4.1, 2.1, 4.1))

# Start new empty plot
plot.new()
plot.window(xlim = c(10000, 40000000), ylim = c(0, 35), log = "x")
box()
axis(1, tck = -0.01, at = c(10000, 100000,1000000,10000000), labels = expression(10^4, 10^5, 10^6, 10^7))
axis(1, tck = -0.01, at = c(20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,20000000,30000000,40000000), labels = F)
y.axis = seq(0, 35, by = 1)
axis(2, las = 1)
axis(2, las = 1, tck = -0.01, at = y.axis, labels = F)

# Load temp data
b = read.table("temp-data.tsv", sep = "\t", header = T)

# Rescale temperature data to match X axis and take the mean of every 2000 year window
temps = zoo(rescale(b$Temp, to = c(0, 35)))
temps.mean = rollapply(temps, width = 2000, FUN = mean, align = "right", by = 200)
times.center = rollapply(zoo(b$Tbp), width = 200, FUN = median, align = "right", by = 200)

# Add temperature points to the plot
points(times.center, temps.mean, col = rgb(1, 0, 0), lwd = 1, type = "l")

# Plot each bootstrapped PSMC file in gray color
for (i in 1:100){
   a = read.table(paste("test.", i, ".txt", sep = ""), header = F)
   points(a$V1 + 0.0001, a$V2, col = "grey", type = "l", lwd = 1)
}

# Plot the overall PSMC results as black line
a = read.table("test.0.txt", header = F)
points(a$V1 + 0.0001, a$V2, col = "black", type = "l", lwd = 3)

# Add lines for each major geomagnetic reversal
abline(v = 780000, lty = "dashed", col = "black")
abline(v = 2580000, lty = "dashed", col = "black")
abline(v = 3580000, lty = "dashed", col = "black")

# Add a second Y axis
x = c(min(b$Temp), seq(from = -15, to = 20, by = 5), max(b$Temp))
axis(4, las = 1, tck = -0.01, at = rescale(x, to = c(0, 35))[2:(length(x) - 1)], labels = seq(from = -15, to = 20, by = 5))

# Write plot to PDF file
dev.off()
```

