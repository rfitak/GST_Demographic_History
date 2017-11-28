# Simulating various demographic histories
This section provides the code for the simulations in our paper using [msHOT-lite](https://github.com/lh3/foreign/tree/master/msHOT-lite).  We simulated three different scenarios, with 10 simulations for each set of parameters.  The parameters were selected to match that from green sea turtles (e.g., genome size), and two chromosomes (1 diploid genome) was sampled at the end of each simulation.  The parameters are further described below. The two simulated scenarios included: 
1.  Two demes (populations) that share migrants duing the mid-Pleistocene transition (MPT)
    - we varied the migration rate in log10 scale from 0-1.
    - MPT occurred from 0.5 - 1.2 million years ago (mya)
2.  Two demes (populations) that share migrants since the Brunhes-Matuyama magnetic polarity reversal
    - the reversal occurred 0.775 mya, and the populations share migrants ever since
3.  A single deme with a mutation rate increase (burst) during the Brunhes-Matuyama magnetic polarity reversal
    - the reversal occurred 0.775 mya, and we included a 10,000 year window flanking the reversal
    - we varied the mutation rate on a log scale from 1% to 1000% the background rate

Here is how the general parameters for the simulations were defined.
First, we generated the ms command to match the observed population history. The program PSMC does this for us:
```bash
# Get ms command from observed PSMC output
psmc2history.pl \
   C_mydas.filtered.combined.psmc | \
   history2ms.pl \
   -n2 \
   -L15000000 \
   -u1.2e-08 \
   -g42.8 \
   -r1
```
Parameters:
- -n2  :: number of chromosomes to simulate (2)
- -L15000000  :: length of each chromosome or sequence (15 megabases)
- -u1.2e-08  :: mutation rate (1.2e-8)
- -g42.8  :: generation time in years (42.8), taken from [Seminoff 2004](http://www.iucnredlist.org/details/4615/0)
- -r1  :: number of replicates (1)

Output:
- Result is an initial Theta (T) = 81649.7075005211
- Result is an initial recombination (r) = 16092.2966989892

## Define simulation parameters:
- Sequence length and number
    - 150 15-megabase fragments per individual = 2.25 GB
    - This is approximately the genome size reported in [Wang et al. 2013 *Nature Genetics*](https://www.nature.com/articles/ng.2615)
    - 2 haploid individuals are sampled at the end, which corresponds to a single, diploid individual's genome
- Extant Effective Population size, N
    - Theta = 81649.7075005211 = 4Nu (definition of theta)
    - N = (81649.7075005211) / (4 \* 1.2e-08 \* 15000000) = 113402
- Historical Effective Population size, No
    - From the ms command produces above, the historical population size is 0.4489 \* N
    - No = 0.4489 \* 113402 = 50,906
- New theta parameter for simulations, (-t), scaled by the new No and by 15 megabase fragment size
    - Theta = 4 \* No \* u
    - Theta = 4 \* 50906 \* 1.2e-08 \* 15000000
    - Theta = 36652.32
- New recombination parameter for simulations, (-r), scaled by the new No and by 15 megabase fragment size
    - r = 4 \* No \* p, where p = recombination rate
    - r from observed output was 16092.2966989892, and must be rescaled to new No
      - 16092.2966989892 = 4 \* 113402 \* p
      - p = 0.03547621889
    - New r = 4 \* 50906 \* 0.03547621889
    - r = 7223.8319878
- Scaling time points
    - all times are scaled by 4No generations
    - MPT = 0.5 - 1.2 mya, or 11,682.2 - 28,037.4 generations in the past (g = 42.8 years)
      - 11,682.2 / (4 \* 50906) = 0.05737 scaled generations
      - 28,037.4 / (4 \* 50906) = 0.13769 scaled generations
    - The Brunhes-Matuyama Reversal occurred 0.775 mya, or 18,107.5 generations in the past (g = 42.8 years)
      - 18,107.5 / (4 \* 50906) = 0.088926 scaled generations
    - The Brunhes-Matuyama Reversal occurred 0.775 mya, and we flank this by 5,000 years on each side for a 10,000 year window from 0.77 - 0.78 mya, or 17,990.7 - 18,224.3 generations in the past (g = 42.8 years)
      - 17,990.7 / (4 \* 50906) = 0.0883523 scaled generations
      - 18,224.3 / (4 \* 50906) = 0.0894997 scaled generations
      
      
## Scenario 1:  Simulating variable migration
We start with the demes (populations) that are independent, and only share migrants during the MPT.  The rate of migration during the MPT between the two demes is varied from 0-1.  In total, 56 migration rates are used, and can be viewed in the file [migrations.list](./Data/migrations.list).  Each simulated migration rate is repeated 10 times, for a total of 560 simulated histories.

```bash
# Setup Folder
mkdir SIM1
cd SIM1/

# Loop through each set of parameters in "migrations.list"
for i in {1..560}
do

# Get the migration rate and replicate number
M=$(sed -n "${i}"p ../migrations.list | cut -d" " -f1)
n=$(sed -n "${i}"p ../migrations.list | cut -d" " -f2)

# Rescale migration rate by 4N
m=$(awk "BEGIN {print ${M}*4*50906}")

# Assign simulation ID
sim="Sim1.${M}.${n}"

# Run ms simulation
msHOT-lite 2 150 \
   -t 36652.32 \
   -r 7223.8319878 \
   15000000 \
   -I 2 2 0 0 \
   -em 0.05737 2 1 $m \
   -l \
   -em 0.13769 2 1 0 > ${sim}.txt
   
# Convert to psmcfa file format
ms2psmcfa.pl ${sim}.txt > ${sim}.psmcfa

# Run PSMC on simulated history
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${sim}.psmc ${sim}.psmcfa

# Generate output plots
psmc_plot.pl \
   -X50000000 \
   -p \
   -g42.8 \
   -R \
   -x10000 \
   -u1.2e-08 \
   ${sim} \
   ${sim}.psmc

# Remove uneeded files
rm -rf $sim.eps $sim.gp $sim.par $sim.pdf $sim.txt

# End loop
done
```


## Scenario 2:  Simulating variable migration, with rates fixed post Brunhes-Matuyama reversal
The scenario is very similar to the first, but rather than migration ONLY during the MPT, this time migration between the two demes is permanent post the Brunhes-Matuyama reversal.  In others, the reversal initiated contact between two populations 775k years ago.  The migrations rates tested are the same.

```bash
# Setup Folder
mkdir SIM2
cd SIM2/

# Loop through each set of parameters in "migrations.list"
for i in {1..560}
do

# Get the migration rate and replicate number
M=$(sed -n "${i}"p ../migrations.list | cut -d" " -f1)
n=$(sed -n "${i}"p ../migrations.list | cut -d" " -f2)

# Rescale migration rate by 4N
m=$(awk "BEGIN {print ${M}*4*50906}")

# Assign simulation ID
sim="Sim2.${M}.${n}"

# Run ms simulation
msHOT-lite 2 150 \
   -t 36652.32 \
   -r 7223.8319878 \
   15000000 \
   -I 2 2 0 $m \
   -l \
   -em 0.088926 2 1 0 > ${sim}.txt
   
# Convert to psmcfa file format
ms2psmcfa.pl ${sim}.txt > ${sim}.psmcfa

# Run PSMC on simulated history
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${sim}.psmc ${sim}.psmcfa

# Generate output plots
psmc_plot.pl \
   -X50000000 \
   -p \
   -g42.8 \
   -R \
   -x10000 \
   -u1.2e-08 \
   ${sim} \
   ${sim}.psmc

# Remove uneeded files
rm -rf $sim.eps $sim.gp $sim.par $sim.pdf $sim.txt

# End loop
done
```


## Scenario 3:  Simulating a mutation rate increase during the Brunhes-Matuyama reversal
In this set of simualtions, we assume a single deme (population) and during the reversal there is more ionizing radiation present and mutation rates are affected.  To encode this, we simply adjust 'theta' during the reversal to account for the increased mutation rate.  The rates are reflected as a percentage of the background rate (1.2e-08). In total, 46 migration rates are used, from 1% to 1000% the background rate.  With 10 replicates per rate, 460 simulations were performed in total.

```bash
# Setup Folder
mkdir SIM3
cd SIM3/

# Loop through each set of parameters in "mutations.list"
for i in {1..460}
do

# Get the mutation rate and replicate number
t=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ../mutations.list | cut -d" " -f1)
n=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ../mutations.list | cut -d" " -f2)

# Assign simulation ID
sim="Sim3.${t}.${n}"

# Run ms simulation
msHOT-lite 2 150 \
   -t 36652.32 \
   -r 7223.8319878 \
   15000000 \
   -l \
   -en 0.0883523 1 $t \
   -en 0.0894997 1 1 > ${sim}.txt
   
# Convert to psmcfa file format
ms2psmcfa.pl ${sim}.txt > ${sim}.psmcfa

# Run PSMC on simulated history
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${sim}.psmc ${sim}.psmcfa

# Generate output plots
psmc_plot.pl \
   -X50000000 \
   -g42.8 \
   -R \
   -x10000 \
   -u1.2e-08 \
   ${sim} \
   ${sim}.psmc

# Remove uneeded files
rm -rf $sim.eps $sim.gp $sim.par $sim.pdf $sim.txt

# End loop
done
```


## Putting it all together
Now we use various commands, especially in R, to calculate the difference between the simulated demographic history and the observed demographic history.  For each simulation, we use the frequency of coalescent events across the 64 atomic time intervals reported by PSMC as a probability distribution. Then, we use a distance metric, D, called the [Jensen-Shannon divergence](http://ieeexplore.ieee.org/document/1207388/).  D ranges from 0 to 1, with 0 being a perfect fit and 1 being completely different.  Simulations with a smaller D are closer in their history to the observed data than larger values of D.

First, we build a table of the probability distributions for each simulation
```bash
# For the observed data
awk '/RD.20/{f=1}/\/\//{f=0}f&&/RS/{print $2,$5,$6}' \
   C_mydas.filtered.psmc | \
   cut -f3 -d" " | \
   tr "\n" "\t" > obs.tMRCAbins.tsv
echo "" >> obs.tMRCAbins.tsv

# For SIM1
cd SIM1
while read i
   do
   i=$(echo "$i" | cut -d" " -f1)
   for n in {1..10}
      do
      a=$(awk '/RD.20/{f=1}/\/\//{f=0}f&&/RS/{print $2,$5,$6}' \
         Sim1.${i}.${n}.psmc | \
         cut -f3 -d" " | \
         tr "\n" "\t")
      echo -e "${i}\t${n}\t${a}" >> ../Sim1.tMRCAbins.tsv
      echo "Finished $i $n"
   done
done < ../migrations.list
cd ..

# For SIM2
cd SIM2
while read i
   do
   i=$(echo "$i" | cut -d" " -f1)
   for n in {1..10}
      do
      a=$(awk '/RD.20/{f=1}/\/\//{f=0}f&&/RS/{print $2,$5,$6}' \
         Sim2.${i}.${n}.psmc | \
         cut -f3 -d" " | \
         tr "\n" "\t")
      echo -e "${i}\t${n}\t${a}" >> ../Sim2.tMRCAbins.tsv
      echo "Finished $i $n"
   done
done < ../migrations.list
cd ..

# For SIM3
cd SIM3
while read i
   do
   i=$(echo "$i" | cut -d" " -f1)
   for n in {1..10}
      do
      a=$(awk '/RD.20/{f=1}/\/\//{f=0}f&&/RS/{print $2,$5,$6}' \
         Sim3.${i}.${n}.psmc | \
         cut -f3 -d" " | \
         tr "\n" "\t")
      echo -e "${i}\t${n}\t${a}" >> ../Sim3.tMRCAbins.tsv
      echo "Finished $i $n"
   done
done < ../mutations.list
cd ..
```

The remaining calculation of D will be performed in R.  First, we need a function that calculates the measure D.  The function below will do this.  It will calculate all pairwise measures of D and output into a matrix.
```R
# JSD
# Function to calculate Jensen-Shannon divergence
   # Arumugam et al. 2011 (Nature)
# "the square root of JSD is a real distance metric (Endres & Schindelin, 2003;
# which is cited as ref. 59 in the supplement of Arumugam et al 2011). 
# This is the reason why we used the square root and not JSD itself.
   # data.dist=dist.JSD(data)

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
	KLD <- function(x,y) sum(x *log(x/y))
	JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
	matrixColSize <- length(colnames(inMatrix))
	matrixRowSize <- length(rownames(inMatrix))
	colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
        
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

	for(i in 1:matrixColSize) {
		for(j in 1:matrixColSize) { 
			resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
			as.vector(inMatrix[,j]))
		}
	}
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
 }
```

Next, we calculate D for each simulation compared with the observed data.

```R
# Load observed distribution
obs = read.table("obs.tMRCAbins.tsv", sep = "\t", header = F)[, 1:64]

# Load simulated distributions (Sim1)
sims = read.table("Sim1.tMRCAbins.tsv", sep = "\t", header = F)[, 1:66]
data = rbind(as.numeric(obs), sims[, 3:66])
dist = as.matrix(dist.JSD(t(data)))
results = dist[, 1]
write(results, file = "Sim1.JS.distances", ncolumns = 1)

# Load simulated distributions (Sim2)
sims = read.table("Sim2.tMRCAbins.tsv", sep = "\t", header = F)[, 1:66]
data = rbind(as.numeric(obs), sims[, 3:66])
dist = as.matrix(dist.JSD(t(data)))
results = dist[, 1]
write(results, file = "Sim2.JS.distances", ncolumns = 1)

# Load simulated distributions (Sim3)
sims = read.table("Sim3.tMRCAbins.tsv", sep = "\t", header = F)[, 1:66]
data = rbind(as.numeric(obs), sims[, 3:66])
dist = as.matrix(dist.JSD(t(data)))
results = dist[, 1]
write(results, file = "Sim3.JS.distances", ncolumns = 1)
```
