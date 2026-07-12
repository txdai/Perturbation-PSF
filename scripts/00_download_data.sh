#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw
curl -L --fail --continue-at - \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274447/suppl/GSE274447_RAW.tar" \
  -o data/raw/GSE274447_RAW.tar

echo "Downloaded GEO GSE274447 raw archive to data/raw/."

