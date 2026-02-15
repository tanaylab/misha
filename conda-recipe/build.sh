#!/bin/bash
set -ex

export DISABLE_AUTOBREW=1

# Remove -shared from Makevars if on macOS (not valid linker flag there)
if [[ "$(uname)" == "Darwin" ]]; then
    sed -i '' 's/-shared//g' src/Makevars
fi

${R} CMD INSTALL --build . ${R_ARGS}
