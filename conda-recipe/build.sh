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

# R 4.0 defaults to C++11 but the codebase requires C++14.  R 4.1+ already
# default to C++14 or higher.  We inject the flag only at conda-build time
# to avoid CRAN policy issues (CXX_STD = CXX14 is deprecated in R-devel).
# Append to CXXFLAGS (not PKG_CXXFLAGS) because GCC uses the last -std= flag
# and CXXFLAGS comes after PKG_CXXFLAGS in R's compile command.
R_MAJOR=$("${R}" --vanilla -s -e 'cat(R.version$major)')
R_MINOR=$("${R}" --vanilla -s -e 'cat(R.version$minor)')
if [[ "${R_MAJOR}" -eq 4 && "${R_MINOR%%.*}" -lt 1 ]]; then
    echo "CXXFLAGS += -std=gnu++14" >> src/Makevars
fi

${R} CMD INSTALL --build . ${R_ARGS}
