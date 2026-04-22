#!/usr/bin/env bash
# fetch_processed_outputs.sh
# downloads processed mosaic outputs listed in the manifest tsv into a local cache.
# skips files already present and rows flagged not_in_public_release.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cache_dir="./cache/mosaic_processed"
manifest="${script_dir}/../assets/processed_outputs_manifest.tsv"
list_only=0

usage() {
  cat <<EOF
usage: $(basename "$0") [--cache_dir <path>] [--manifest <path>] [--list-only]
  --cache_dir   where to drop files (default: ./cache/mosaic_processed)
  --manifest    tsv manifest (default: pipeline/assets/processed_outputs_manifest.tsv)
  --list-only   print manifest entries, don't download
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cache_dir) cache_dir="$2"; shift 2 ;;
    --manifest)  manifest="$2"; shift 2 ;;
    --list-only) list_only=1; shift ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -f "$manifest" ]] || { echo "manifest not found: $manifest" >&2; exit 1; }
mkdir -p "$cache_dir"

fetched=0; skipped=0; failed=0; flagged=0; bundled=0

# awk reformats each row with a sentinel delimiter and escapes empty fields,
# since bash read collapses adjacent tabs (tab is IFS whitespace).
while IFS='|' read -r filename url required_by notes; do
  [[ -z "${filename:-}" || "$filename" == "filename" ]] && continue
  # unescape sentinel placeholder
  [[ "$url" == "__EMPTY__" ]] && url=""
  [[ "$notes" == "__EMPTY__" ]] && notes=""
  [[ "$required_by" == "__EMPTY__" ]] && required_by=""

  if [[ "$list_only" -eq 1 ]]; then
    printf "  %-60s  required_by=%s  notes=%s\n" "$filename" "${required_by:-}" "${notes:-}"
    continue
  fi

  if [[ "${notes:-}" == local_in_repo_* ]]; then
    echo "bundled in repo: $filename (notes=${notes})"
    bundled=$((bundled+1)); continue
  fi

  if [[ "${notes:-}" == *"not_in_public_release"* ]]; then
    echo "warning: $filename flagged not_in_public_release — skipping"
    flagged=$((flagged+1)); continue
  fi

  dest="${cache_dir}/${filename}"
  if [[ -f "$dest" ]]; then
    echo "cached: $filename"
    skipped=$((skipped+1)); continue
  fi

  if [[ -z "${url:-}" ]]; then
    echo "warning: $filename has no url (notes=${notes:-}) — skipping"
    flagged=$((flagged+1)); continue
  fi

  echo "fetching: $filename"
  if curl -L --fail --silent --show-error -o "$dest" "$url"; then
    fetched=$((fetched+1))
  else
    echo "failure: $filename" || true
    rm -f "$dest"
    failed=$((failed+1))
  fi
done < <(awk -F'\t' 'BEGIN{OFS="|"} {
  for (i=1; i<=4; i++) if ($i == "") $i = "__EMPTY__";
  print $1, $2, $3, $4
}' "$manifest")

if [[ "$list_only" -eq 0 ]]; then
  echo "summary: ${fetched} fetched, ${skipped} cached/skipped, ${failed} failed, ${flagged} flagged_not_public, ${bundled} bundled_local"
fi
