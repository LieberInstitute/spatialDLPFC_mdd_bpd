#!/usr/bin/env bash
set -euo pipefail

## run the published TOP-LD API client against a predeclared pair list.
base=${1:?usage: 19_run_topld_validation.sh BASE}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
query="$script_dir/susie-topld-pairs.txt"
out="$base/topld_validation"
mkdir -p "$out"

## pin the reviewed client artifact and reject any upstream replacement.
client_url=https://raw.githubusercontent.com/linnabrown/topld_api/main/topld_api
expected_sha=8442e10383b424c34f6d00930a99bdf2b83ecd340e940d17a1f162a6180eeea5
client=$(mktemp)
trap 'rm -f "$client"' EXIT
curl -fsSL "$client_url" -o "$client"
observed_sha=$(sha256sum "$client" | cut -d ' ' -f1)
if [[ "$observed_sha" != "$expected_sha" ]]; then
    echo "TOP-LD client checksum changed" >&2
    exit 1
fi
chmod 700 "$client"

## EUR matches the primary ancestry; zero threshold requests the specified pairs.
timeout 180 "$client" -thres 0.0 -pop EUR -maf 0.001 -inFile "$query" \
    -outputLD "$out/topld_EUR_LD.tsv" -outputInfo "$out/topld_EUR_info.tsv"

## distinguish returned LD estimates from pairs absent in the service response.
Rscript - "$query" "$out" "$client_url" "$observed_sha" <<'RS'
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
query_file <- args[[1L]]
out_dir <- args[[2L]]
requests <- fread(query_file, header = FALSE, sep = ",", col.names = c("marker1", "marker2"))
info <- fread(file.path(out_dir, "topld_EUR_info.tsv"))
ld <- fread(file.path(out_dir, "topld_EUR_LD.tsv"))

normalize_marker <- function(marker) {
  if (!startsWith(marker, "chr")) return(marker)
  fields <- strsplit(marker, ":", fixed = TRUE)[[1L]]
  pos <- as.integer(fields[[2L]])
  ref <- fields[[3L]]
  alt <- fields[[4L]]
  hit <- info[Position == pos & REF == ref & ALT == alt, unique(rsID)]
  if (length(hit) == 1L) hit else marker
}
requests[, `:=`(
  normalized1 = vapply(marker1, normalize_marker, character(1L)),
  normalized2 = vapply(marker2, normalize_marker, character(1L))
)]
requests[, pair_key := vapply(seq_len(.N), function(i) {
  paste(sort(c(normalized1[[i]], normalized2[[i]])), collapse = "|")
}, character(1L))]
ld[, pair_key := vapply(seq_len(.N), function(i) {
  paste(sort(c(rsID1[[i]], rsID2[[i]])), collapse = "|")
}, character(1L))]
audit <- merge(requests, ld, by = "pair_key", all.x = TRUE)
audit[, status := fifelse(is.na(R2), "not_returned_by_TOP-LD", "returned")]
fwrite(audit, file.path(out_dir, "topld_pair_audit.tsv"), sep = "\t")

provenance <- data.table(
  field = c("client_url", "client_sha256", "population", "r2_threshold", "maf_threshold",
            "requested_pairs", "returned_pairs", "run_time"),
  value = c(args[[3L]], args[[4L]], "EUR", "0.0", "0.001", nrow(requests),
            nrow(ld), format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
)
fwrite(provenance, file.path(out_dir, "topld_provenance.tsv"), sep = "\t")
RS
