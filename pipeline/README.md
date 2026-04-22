# mosaic pipeline

mosaic calls microsatellite instability (msi) status from paired tumor/normal exomes. this pipeline wraps the msings-based preprocessing and the trained classifier methods from the paper.

## quickstart: one paired tumor/normal exome

```bash
# 1. clone and cd
git clone https://github.com/ronaldhause/mosaic
cd mosaic/pipeline

# 2. build container
docker build -t mosaic:latest containers/

# 3. create a one-row samplesheet
cat > my_sample.csv <<EOF
sample_name,tumor_bam,normal_bam,tumor_type
my_patient_01,/data/tumor.bam,/data/normal.bam,COAD
EOF

# 4. run with nextflow (local, docker)
nextflow run main.nf \
  -profile local,docker \
  --samplesheet my_sample.csv \
  --reference /refs/hg19.fa \
  --microsatellite_bed /refs/microsatellites.bed \
  --outdir ./results

# 5. check the result
cat results/mosaic_results.csv
```

make equivalent:

```bash
make all SAMPLESHEET=my_sample.csv REFERENCE=/refs/hg19.fa MSBED=/refs/microsatellites.bed
```

for slurm:

```bash
nextflow run main.nf -profile slurm,singularity \
  --samplesheet my_sample.csv \
  --reference /refs/hg19.fa \
  --microsatellite_bed /refs/microsatellites.bed \
  --outdir ./results
```

## prerequisites

- aligned paired bams (sorted, indexed, hg19)
- reference fasta with matching `.fai`
- microsatellite bed — the 516,876 loci used in the paper, available from the krishna lab url, or regenerate with `build_reference.nf` (see below). if you regenerate, set `--num_loci` to match the number of loci in your bed.
- msings (installed locally or via the container image)
- rough resource expectations for one paired sample (~50x exome): ~2 gb msings output, ~1-4 cpu-hours
- hg38 note: the model was trained on hg19 coordinates. for hg38 samples you need to liftover the bed, and the defb locus parameter changes — verify the new locus id matches an entry in your bed before running.

## building the microsatellite bed from scratch

if you don't have the published 516,876-locus bed, `build_reference.nf` regenerates one from an hg19 fasta and the bait bed of whatever capture kit produced your exomes. it runs the full method from the paper: misa → +/- 5 bp pad → merge adjacent (<=10 bp gap) into compound/complex loci → restrict to loci within 50 bp of capture baits → optional annovar annotation.

```bash
# one-time: build the reference bed.
nextflow run pipeline/build_reference.nf \
  -profile local,docker \
  --reference /refs/hg19.fa \
  --capture_bait_bed /refs/nimblegen_seqcap_ez_exome_v3.bed \
  --outdir ./ref_build

# then run mosaic against it.
nextflow run pipeline/main.nf \
  -profile local,docker \
  --samplesheet my_sample.csv \
  --reference /refs/hg19.fa \
  --microsatellite_bed ./ref_build/microsatellites.bed \
  --num_loci $(wc -l < ./ref_build/microsatellites.bed) \
  --outdir ./results
```

notes:

- the bait bed has to match the kit that captured your exomes. the paper used nimblegen seqcap_ez_exome_v3; anything else will produce a different locus count but still a self-consistent reference for your samples.
- annovar is optional. without it, `locus_annotation.tsv` still carries the structural fields (repeat_type, subunits, is_compound, is_complex, etc.) — enough for downstream use. mosaic prediction itself doesn't depend on annovar annotations.
- if you have annovar installed, pass `--annovar_dir /path/to/annovar --annovar_db /path/to/annovar/humandb`. the container does not include annovar because it requires personal registration.
- full parameter reference in [assets/build_reference_params.md](assets/build_reference_params.md).

## running mosaic on new samples

samplesheet columns:

- `sample_name` — unique id per pair
- `tumor_bam` — path to tumor bam (indexed)
- `normal_bam` — path to matched normal bam (indexed)
- `tumor_type` — optional. annotation only; not used by the classifier.

output columns in `mosaic_results.csv`:

```
sample_name, tumor_type, peak_avg, peak_sd, num_unstable, num_called, prop_unstable, defb_status, msi_status
```

how to interpret:

- `msi_status` — MSI-H or MSS. NA if `num_called == 0`.
- `peak_avg` — mean gain in alleles (tumor vs normal) across all called loci. thresholds:
  - `peak_avg >= 0.0055` → MSI-H
  - `0.0029 < peak_avg < 0.0055` and defb locus unstable → MSI-H
  - otherwise MSS
- `prop_unstable` — proportion of loci called unstable. a useful continuous burden metric even for MSS-called samples.
- `defb_status` — status of the defb105a/b locus (chr8:7679723-7679741), used as a tiebreaker in the borderline band.

common failure modes:

- `num_called == 0` — tumor or normal bam has no coverage over the microsatellite loci. check bam coverage and bed contig naming (chr1 vs 1).
- missing defb locus — defaults to "stable". if many samples show this, your bed is probably missing the defb entry.

adjusting thresholds: `--peak_avg_threshold_high` and `--peak_avg_threshold_low` via nextflow params, or `THRESHOLD_HIGH=` / `THRESHOLD_LOW=` for make.

## the bundled classifier

the trained classifier from the paper is bundled at [assets/mosaic_classifier_063016.robj](assets/mosaic_classifier_063016.robj) (a caret `train` object wrapping rpart). the pipeline's `classify_msi.R` replicates it inline as a two-threshold decision tree — no runtime load needed. if you want to call the model directly instead:

```r
load("pipeline/assets/mosaic_classifier_063016.robj")  # provides `mosaic`
predict(mosaic, <your_data_frame>, type="raw")
```

input columns expected by the model: `peak_avg` and `X8.7679723.7679741` (defb locus binary, 1 = unstable).
