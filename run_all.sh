#!/usr/bin/env bash
#
# run_all.sh — reproduce every analysis and figure for
#   "Are the phases real? Distinguishing bounded interaction groups from
#    spatially structured drift in central Mississippi Valley decorated ceramics"
#
# Usage:
#   ./run_all.sh                 # uses `python3`
#   PYTHON=.venv/bin/python ./run_all.sh
#
# Prerequisites (see README.md):
#   pip install -r requirements.txt
#   pip install -e .             # makes the `mls_emergence` package importable
#
# The numbered scripts run in ascending order; outputs land in output/ and
# figures/. The run continues past any single failure and reports the count at
# the end. Two scripts (08, 10) also need the optional monument-mls package and
# will exit with a clear message if it is absent; this is expected and does not
# affect any main-text result.

set -u
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"

mkdir -p output figures
echo "Interpreter: $($PY --version 2>&1)"
echo

fail=0
run() {
  echo "=================================================================="
  echo ">>> $1"
  if ! "$PY" "$1"; then
    echo "!!! $1 exited non-zero (see message above)"
    fail=$((fail + 1))
  fi
}

# Numbered analysis + figure scripts, in order.
for s in analyses/[0-9][0-9]_*.py; do
  run "$s"
done

# House-style figure pipeline (most main-text and supplemental figures).
run "analyses/make_figures.py"

echo "=================================================================="
echo "Done. ${fail} script(s) reported a non-zero exit."
echo "(Scripts 08 and 10 exit non-zero without the optional monument-mls"
echo " package; that is expected and does not affect any main-text result.)"
