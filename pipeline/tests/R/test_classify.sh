#!/usr/bin/env bash
# end-to-end test of classify_msi.R against a fixture features csv.
# uses the three decision-tree branches: msi-h high, msi-h borderline+defb,
# mss, and na (num_called==0).

set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bin="$here/../../bin"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

mkdir -p "$tmp/features"
cat > "$tmp/features/s1.features.csv" <<'EOF'
sample_name,tumor_type,peak_avg,peak_sd,num_unstable,num_called,prop_unstable,defb_status
msi_high,COAD,0.01,0.02,5000,500000,0.01,unstable
borderline_defb_unstable,UCEC,0.004,0.01,2000,500000,0.004,unstable
borderline_defb_stable,STAD,0.004,0.01,2000,500000,0.004,stable
mss,BRCA,0.001,0.005,500,500000,0.001,stable
no_coverage,LUAD,NA,NA,0,0,NA,stable
EOF

Rscript "$bin/classify_msi.R" \
    --input_dir "$tmp/features" \
    --pattern "*.features.csv" \
    --output "$tmp/results.csv"

# sanity: every expected status is present
awk -F, 'NR>1 {gsub(/"/,""); print $1","$NF}' "$tmp/results.csv" > "$tmp/status.csv"

expect() {
    grep -qx "$1" "$tmp/status.csv" || { echo "FAIL: missing row '$1'"; cat "$tmp/status.csv"; exit 1; }
}

expect "msi_high,MSI-H"
expect "borderline_defb_unstable,MSI-H"
expect "borderline_defb_stable,MSS"
expect "mss,MSS"
expect "no_coverage,NA"

echo "classify_msi.R: PASS"
