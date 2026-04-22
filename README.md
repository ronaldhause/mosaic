### MicrOSAtellite Instability Classifier (MOSAIC)

#### Summary
Scripts and pipeline for calling microsatellite instability (MSI) status from paired tumor/normal exomes, accompanying "Classification and characterization of microsatellite instability across 18 cancer types" in [Nature Medicine](http://www.nature.com/nm/index.html).

## Abstract

Microsatellite instability (MSI), the spontaneous loss or gain of nucleotides from repetitive DNA tracts, is a diagnostic phenotype for gastrointestinal, endometrial, and colorectal tumors, yet the landscape of instability events across a wider variety of cancer types remains poorly understood. To explore MSI across malignancies, we examined 5,930 cancer exomes from 18 cancer types at more than 200,000 microsatellite loci and constructed a genomic classifier for MSI. We identified MSI-positive tumors in 14 of the 18 cancer types. We also identified loci that were more likely to be unstable in particular cancer types, resulting in specific instability signatures that involved cancer-associated genes, suggesting that instability patterns reflect selective pressures and can potentially identify novel cancer drivers. We also observed a correlation between survival outcomes and the overall burden of unstable microsatellites, suggesting that MSI may be a continuous, rather than discrete, phenotype that is informative across cancer types. These analyses offer insight into conserved and cancer-specific properties of MSI and reveal opportunities for improved methods of clinical MSI diagnosis and cancer gene discovery.

## Install
```bash
git clone https://github.com/ronaldhause/mosaic
cd mosaic

# Option 1: Docker (recommended)
docker build -t mosaic:latest pipeline/containers/

# Option 2: Conda
conda env create -f pipeline/containers/environment.yml
conda activate mosaic
```

Prerequisites: Docker or Conda, Nextflow (for the end-to-end pipeline), and paired tumor/normal exome BAMs aligned to hg19.

## Running MOSAIC on Your Own Data

### End-to-end pipeline

The preprocessing and classification pipeline is in [pipeline/](pipeline/); see [pipeline/README.md](pipeline/README.md) for a full walkthrough. A minimal Nextflow invocation:

```bash
cd pipeline
nextflow run main.nf \
  -profile local,docker \
  --samplesheet my_sample.csv \
  --reference /refs/hg19.fa \
  --microsatellite_bed /refs/microsatellites.bed \
  --outdir ./results
```

The samplesheet lists paired tumor/normal BAMs. The pipeline performs mSINGS-based preprocessing at each microsatellite locus and then applies the trained MOSAIC classifier to produce MSI-H/MSS calls.

### Using MOSAIC as a classifier only

If unstable microsatellite calls have already been generated with a different tool (for example, [mSINGS](https://bitbucket.org/uwlabmed/msings), [lobSTR](http://lobstr.teamerlich.org/), or [HipSTR](https://hipstr-tool.github.io/HipSTR/)), the preprocessing stage can be skipped and the trained MOSAIC classifier applied directly in R:

1. Call microsatellite markers in hg18 or hg19 using [MISA](https://webblast.ipk-gatersleben.de/misa/), or download them from the [UCSC Genome Browser](http://genome.ucsc.edu/).
2. Run mSINGS, lobSTR, HipSTR, or an equivalent tool to identify unstable microsatellites from paired tumor/normal exome BAMs at each microsatellite site.
3. Reformat the output into a table containing, per sample: an identifier for each locus (`msi`), the average gain in unique alleles in tumor relative to matched normal tissue across all interrogated microsatellites (`peak_avg`), and a binary variable indicating whether the DEFB105A/B locus at chr. 8:7679723–7679741 (`X8.7679723.7679741`) is unstable.
4. Predict MSI classes using the bundled classifier:

```r
load("pipeline/assets/mosaic_classifier_063016.robj")
predict(mosaic, your_data, type = "raw")
```

## Data

Primary TCGA exome alignments are controlled-access and available from the [TCGA Research Network](http://cancergenomenih.gov/). Primary and processed MSI calls derived from these alignments, the MOSAIC classifier itself, and intermediate results and tables are available [here](http://krishna.gs.washington.edu/content/members/hauser/mosaic/).

## Dependencies

Perl, Python, and R were used across the primary and secondary analyses, including the ggplot2, rpart, clusterProfiler, qvalue, and survival libraries, among others. External tools include MISA, mSINGS, and Grid Engine.

## Contact

Please let Ron Hause <ronaldhause@gmail.com> or Steve Salipante <stevesal@uw.edu> know if anything is missing or if you have any questions. Apologies in advance for any chaos.
