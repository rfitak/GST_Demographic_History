# Simulating various demographic histories
This section provides the code for the simulations in our paper using [msHOT-lite](https://github.com/lh3/foreign/tree/master/msHOT-lite).  We simulated three different scenarios, with 10 simulations for each set of parameters.  The parameters were selected to match that from green sea turtles (e.g., genome size), and two chromosomes (1 diploid genome) was sampled at the end of each simulation.  The parameters are further described below. The two simulated scenarios included: 
1.  Two demes (populations) that share migrants duing the mid-Pleistocene transition (MPT)
    - we varied the migration rate in log10 scale from 0-1.
    - MPT occurred from 0.5 - 1.2 million years ago (mya)
2.  Two demes (populations) that share migrants since the Brunhes-Matuyama magnetic polarity reversal
    - the reversal occurred 0.78 mya, and the populations share migrants ever since
3.  A single deme with a mutation rate increase (burst) during the Brunhes-Matuyama magnetic polarity reversal
    - the reversal occurred 0.78 mya, and we included a 10,000 year window flanking the reversal
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
- New theta parameter for simulations, (-t), scaled by the new No
    - Theta = 4 \* No \* u




## Scenario 1:  Simulating variable migration


```bash
cd /work/frr6/SIMS/

M=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ../migrations.list | cut -d" " -f1)
n=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ../migrations.list | cut -d" " -f2)

m=$(awk "BEGIN {print ${M}*4*50906}")
sim="Sim1.${M}.${n}"

msHOT-lite 2 150 \
   -t 36652.32 \
   -r 7223.8319878 \
   15000000 \
   -I 2 2 0 0 \
   -em 0.05737 2 1 $m \
   -l \
   -em 0.13769 2 1 0 > ${sim}.txt
   
ms2psmcfa.pl ${sim}.txt > ${sim}.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${sim}.psmc ${sim}.psmcfa

psmc_plot.pl \
   -X50000000 \
   -p \
   -g42.8 \
   -R \
   -x10000 \
   -u1.2e-08 \
   ${sim} \
   ${sim}.psmc

rm -rf $sim.eps $sim.gp $sim.par $sim.pdf $sim.txt

```
