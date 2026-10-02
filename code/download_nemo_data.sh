#!/usr/bin/env bash
set -euo pipefail

data_root=${DATA_ROOT:-/scratch/nemo_download}

download_rdata() {
  local destination=$1
  local source_url=$2

  mkdir -p "$destination"
  wget \
    --recursive \
    --no-parent \
    --no-directories \
    --accept='*.rda,*.Rdata' \
    --directory-prefix="$destination" \
    "$source_url"
}

base_url=https://data.nemoarchive.org/bican/grant/BICAN_Neuromodulation/aibs/transcriptome/nuclei

download_rdata "$data_root/10x" \
  "$base_url/10x_v4/house_mouse/processed/counts/"
download_rdata "$data_root/10x" \
  "$base_url/10x_v4/house_mouse/processed/cell_metrics/"
download_rdata "$data_root/SS" \
  "$base_url/SSv4/house_mouse/processed/counts/"
download_rdata "$data_root/SS" \
  "$base_url/SSv4/house_mouse/processed/cell_metrics/"
