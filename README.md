# Green Sea Turtle Demographic History
The code here is to reconstruct the green sea turtle's demographic history from the genome sequence. The genome sequence was first described by [Wang et al. 2013 *Nature Genetics*](https://www.nature.com/articles/ng.2615) and the raw sequencing data for the genome can be accessed at [PRJNA104937](https://www.ncbi.nlm.nih.gov/bioproject/104937). Furthermore, the code for a series of simulated histories under various migration scenarios is made available. Please follow the links below to get the code and associated descriptions for various analyses and data processing. The results of this study have been published in:  

Fitak, R & Johnsen, S (in press): Green sea turtle (*Chelonia mydas*) population history indicates important demographic changes near the mid-Pleistocene transition. *Marine Biology*, doi: [10.1007/s00227-018-3366-3](https://doi.org/10.1007/s00227-018-3366-3)  

Most of the data, including VCF files of SNPs and results of the simulations are available at:  

Fitak, R & Johnsen, S (2018): Single nucleotide polymorphisms and PSMC results in the green sea turtle *Chelonia mydas*. PANGAEA, [https://doi.pangaea.de/10.1594/PANGAEA.890759](https://doi.pangaea.de/10.1594/PANGAEA.890759)  

Some ancillary data files and custom scripts are also available in the [Data](./Data) folder.

1.  [Download and process raw sequencing data](./GST-data-processing.md)
    - Downloading, trimming, and mapping raw sequencing data.
2. [SNP calling and PSMC analysis](./SNPs-and-PSMC.md)
    - Identifcation and filtering of SNPs, then demographic reconstruction using PSMC
3. [Simulating demographic history](./sims.md)
    - SImulating the various demographic scenarios and comparing with the observed demographic history

## All code and content herein is licensed under:
## [GNU General Public License v2.0](./LICENSE)
