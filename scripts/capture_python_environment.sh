#!/usr/bin/env bash
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python -m pip freeze > "$HERE/provenance/python_packages_at_run.txt"
echo "Python env captured; review before public upload (may include private package URLs)."
