# mosaic pipeline — make-based alternative to nextflow.
# runs the same three steps: mSINGS → compute_features → classify_msi.
# use `make -jN` for parallel mSINGS execution.

SAMPLESHEET ?= samples.csv
REFERENCE   ?= hg19.fa
MSBED       ?= microsatellites.bed
OUTDIR      ?= results
NUM_LOCI    ?= 516876
DEFB_LOCUS  ?= 8:7679723-7679741

BIN_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))bin

# validated samplesheet (tumor_type filled in, bams checked)
VALIDATED := $(OUTDIR)/samplesheet.validated.csv

# derive sample names from column 1 (skip header). requires validated samplesheet.
SAMPLES := $(shell tail -n +2 $(SAMPLESHEET) 2>/dev/null | cut -d, -f1)

MSINGS_TSVS   := $(addprefix $(OUTDIR)/msings/,$(addsuffix .msings.tsv,$(SAMPLES)))
FEATURE_CSVS  := $(addprefix $(OUTDIR)/features/,$(addsuffix .features.csv,$(SAMPLES)))
RESULTS       := $(OUTDIR)/mosaic_results.csv

.PHONY: all validate clean
.SECONDARY:

all: $(RESULTS)

validate: $(VALIDATED)

$(VALIDATED): $(SAMPLESHEET) | $(OUTDIR)
	python3 $(BIN_DIR)/validate_samplesheet.py $< > $@

# look up a column for a sample from the samplesheet.
# $(1) = sample, $(2) = column index (1-based)
sample_col = $(shell awk -F, -v s="$(1)" -v c=$(2) 'NR>1 && $$1==s {print $$c; exit}' $(SAMPLESHEET))

# step 1: run msings per sample. wraps the mSINGS CLI — adjust invocation
# to match your installed version. tumor/normal bams are looked up from
# the samplesheet by sample name.
$(OUTDIR)/msings/%.msings.tsv: $(SAMPLESHEET) | $(OUTDIR)/msings
	@tbam=$$(awk -F, -v s="$*" 'NR>1 && $$1==s {print $$2; exit}' $<); \
	nbam=$$(awk -F, -v s="$*" 'NR>1 && $$1==s {print $$3; exit}' $<); \
	if [ -z "$$tbam" ] || [ -z "$$nbam" ]; then \
	  echo "sample $* not found in $<" >&2; exit 1; \
	fi; \
	echo "[msings] $* tumor=$$tbam normal=$$nbam"; \
	msi analyzer --bam "$$tbam" --normal-bam "$$nbam" \
	    --bed $(MSBED) --fasta $(REFERENCE) --out $@

# step 2: compute per-sample features
$(OUTDIR)/features/%.features.csv: $(OUTDIR)/msings/%.msings.tsv $(SAMPLESHEET) | $(OUTDIR)/features
	@ttype=$$(awk -F, -v s="$*" 'NR>1 && $$1==s {print ($$4==""?"UNKNOWN":$$4); exit}' $(SAMPLESHEET)); \
	Rscript $(BIN_DIR)/compute_features.R \
	    --input $< \
	    --sample_name $* \
	    --tumor_type $${ttype:-UNKNOWN} \
	    --num_loci $(NUM_LOCI) \
	    --defb_locus "$(DEFB_LOCUS)" \
	    --output $@

# step 3: classify — gather step over all feature files
$(RESULTS): $(FEATURE_CSVS) | $(OUTDIR)
	Rscript $(BIN_DIR)/classify_msi.R \
	    --input_dir $(OUTDIR)/features \
	    --pattern "*.features.csv" \
	    --threshold_high 0.0055 \
	    --threshold_low 0.0029 \
	    --output $@

$(OUTDIR) $(OUTDIR)/msings $(OUTDIR)/features:
	mkdir -p $@

clean:
	rm -rf $(OUTDIR)
