# build_reference.nf parameters

rebuilds the microsatellite bed consumed by `main.nf` from an hg19 fasta
and a capture-bait bed. follows the paper's method (hause et al.,
nature medicine 2016).

## required

| param | description |
|---|---|
| `--reference`        | hg19 fasta. a matching `.fai` must sit next to it (run `samtools faidx` if missing). |
| `--capture_bait_bed` | bait bed for the capture kit used to generate the exomes. anything non-nimblegen is fine; the pipeline just needs regions where coverage is expected. |

## optional (with defaults)

| param | default | description |
|---|---|---|
| `--outdir`          | `./ref_build` | publish dir for microsatellites.bed, locus_annotation.tsv, build_summary.tsv. |
| `--misa_motifs`     | `1-5 2-5 3-5 4-5 5-5` | misa `definition(unit_size,min_repeats)` spec. means subunit lengths 1-5 bp, min 5 repeats. |
| `--misa_interrupts` | `0`  | misa interruption param. 0 = perfect repeats only; paper's downstream merge step handles compound/complex assembly. |
| `--pad_bp`          | `5`  | bp added to each side of every raw misa tract before merging. |
| `--merge_gap_bp`    | `10` | max gap (bp) between two adjacent tracts for them to merge into a single compound (same subunit length) or complex (`c*`, different subunit lengths) locus. |
| `--bait_slop_bp`    | `50` | bp slop added to each side of every bait interval before intersecting with the merged loci. preserves loci that sit just off the bait edge. |

## optional — annovar (off by default)

| param | default | description |
|---|---|---|
| `--annovar_dir`   | `null` | path to an annovar install containing `annotate_variation.pl`. |
| `--annovar_db`    | `null` | annovar db dir (e.g. `/opt/annovar/humandb`). refgene must be present. |
| `--annovar_build` | `hg19` | annovar genome build. |

annovar requires personal registration and is not bundled in the container. without it, `locus_annotation.tsv` still includes the structural fields (`repeat_type`, `repeat_subunits`, `subunit_lengths`, `n_repeats`, `is_compound`, `is_complex`, constituent `members`) derived from misa + the merge step. mosaic prediction does not depend on annovar annotations — those are for downstream inspection.

## outputs

| file | description |
|---|---|
| `microsatellites.bed` | 6-col bed (chrom, start, end, locus_id, score, strand). feed this to `main.nf` via `--microsatellite_bed`. |
| `locus_annotation.tsv`| per-locus attributes: repeat_type, repeat_subunits, subunit_lengths, n_repeats, is_compound, is_complex, members, plus genomic_class/gene when annovar is wired in. |
| `build_summary.tsv`   | counts: raw vs merged loci, compound/complex breakdown, defb locus presence check. |

## chained usage

```bash
# one-time: build the reference.
nextflow run pipeline/build_reference.nf \
  -profile local,docker \
  --reference /refs/hg19.fa \
  --capture_bait_bed /refs/nimblegen_seqcap_ez_exome_v3.bed \
  --outdir ./ref_build

# then call msi on new samples.
nextflow run pipeline/main.nf \
  -profile local,docker \
  --samplesheet my_samples.csv \
  --reference /refs/hg19.fa \
  --microsatellite_bed ./ref_build/microsatellites.bed \
  --num_loci $(wc -l < ./ref_build/microsatellites.bed) \
  --outdir ./results
```
