#!/bin/bash

set -euo pipefail

date
echo
echo "Running workflow '${TOOL_REPO}/modules/dada2.nf' from ${PWD}"
echo

echo "Workflow Parameters:"
cat ._wb/tool/params.json
echo

# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${TOOL_REPO}/modules/dada2.nf" \
    -params-file ._wb/tool/params.json \
    -with-report \
    -with-trace \
    -resume

echo
date
echo Done