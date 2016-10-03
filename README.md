### MicrOSAtellite Instability Classifier (MOSAIC)

#### Summary
This page contains the scripts used for secondary analysis and plotting in "Classification and characterization of microsatellite instability across 18 cancer types" found in [Nature Medicine](http://www.nature.com/nm/index.html).

## Abstract

Microsatellite instability (MSI), the spontaneous loss or gain of nucleotides from repetitive DNA tracts, is a diagnostic phenotype for gastrointestinal, endometrial, and colorectal tumors, yet the landscape of instability events across a wider variety of cancer types remains poorly understood. To explore MSI across malignancies, we examined 5,930 cancer exomes from 18 cancer types at more than 200,000 microsatellite loci and constructed a genomic classifier for MSI. We identified MSI-positive tumors in 14 of the 18 cancer types. We also identified loci that were more likely to be unstable in particular cancer types, resulting in specific instability signatures that involved cancer-associated genes, suggesting that instability patterns reflect selective pressures and can potentially identify novel cancer drivers. We also observed a correlation between survival outcomes and the overall burden of unstable microsatellites, suggesting that MSI may be a continuous, rather than discrete, phenotype that is informative across cancer types. These analyses offer insight into conserved and cancer-specific properties of MSI and reveal opportunities for improved methods of clinical MSI diagnosis and cancer gene discovery.

## Data
Primary exome sequence alignments can be downloaded from the [TCGA Research Network](http://cancergenomenih.gov/). Primary and processed MSI calls derived from these exome alignments, the MOSAIC classifier itself, and intermediate results and tables are available [here](http://krishna.gs.washington.edu/content/members/hauser/mosaic/).

## Dependencies
Perl, Python, and R languages were used for various steps throughout the primary and secondary analyses, including the ggplot2, rpart, clusterProfiler, qvalue, and survival libraries, among others. We also used several external computational tools, such as MISA, mSINGS, and Grid Engine.

## Contact
Please let Ron Hause <ronaldhause@gmail.com> or Steve Salipante <stevesal@uw.edu> know if anything is missing or if you have any questions. Our apologies ahead of time for any and all chaos. :)