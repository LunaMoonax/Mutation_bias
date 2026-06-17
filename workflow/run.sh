#!/usr/bin/env bash

set -euo pipefail

snakemake all -p --use-conda --cores 16 --rerun-incomplete