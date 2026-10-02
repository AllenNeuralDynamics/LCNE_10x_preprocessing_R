#!/usr/bin/env bash
set -euo pipefail

data_root=${DATA_ROOT:-/data}

download_rdata() {
  local destination=$1
  local source_url=$2

  mkdir -p "$destination"
  wget \
    --recursive \
    --no-parent \
    --no-directories \
    --accept='*.rda,*.rdata' \
    --directory-prefix="$destination" \
    "$source_url"
}

base_url=https://data.nemoarchive.org/bican/grant/BICAN_Neuromodulation/aibs/transcriptome/nuclei

download_rdata "$data_root/LCv2" \
  "$base_url/10x_v4/house_mouse/processed/counts/"
download_rdata "$data_root/LCv2" \
  "$base_url/10x_v4/house_mouse/processed/cell_metrics/"
download_rdata "$data_root/LCNE_smartseq_raw_NeMO" \
  "$base_url/SSv4/house_mouse/processed/counts/"
download_rdata "$data_root/LCNE_smartseq_raw_NeMO" \
  "$base_url/SSv4/house_mouse/processed/cell_metrics/"
