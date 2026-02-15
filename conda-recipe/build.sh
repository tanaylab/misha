#!/bin/bash
set -ex

export DISABLE_AUTOBREW=1

# Remove -shared from Makevars — it's redundant (R builds shared libs by
# default) and causes linker failures in conda's sysroot on Linux and macOS.
if [[ "$(uname)" == "Darwin" ]]; then
    sed -i '' 's/-shared//g' src/Makevars
else
    sed -i 's/-shared//g' src/Makevars
fi

${R} CMD INSTALL --build . ${R_ARGS}
