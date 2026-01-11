set -euo pipefail

snakemake --use-conda --conda-frontend mamba --cores 80 -p

