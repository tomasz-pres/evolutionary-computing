#!/bin/sh
# Usage: ./run.sh [DATA_DIR]   (default ../lab1)
set -e
cd "$(dirname "$0")"
DATA="${1:-../lab1}"
g++ -O2 -std=c++17 -o ls ls.cpp
./ls --data-dir "$DATA" --out results
python3 report.py --data-dir "$DATA"
