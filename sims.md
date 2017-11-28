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

# Loop through each set of parameters in "migrations.list"
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
