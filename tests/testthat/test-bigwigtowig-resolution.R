# The bundled bigWigToWig is a Linux x86-64 binary, so bigWig import cannot work
# on other platforms without an externally installed converter. There is no
# bigWig import test, which is why this went unnoticed on the macOS CI runner.

test_that("options(misha.bigWigToWig) wins over everything else", {
    withr::local_options(misha.bigWigToWig = "/some/where/bigWigToWig")
    expect_equal(get_bigWigToWig_bin(), "/some/where/bigWigToWig")
    # honored on every platform, including the bundled-binary one
    expect_equal(get_bigWigToWig_bin("Linux"), "/some/where/bigWigToWig")
})

test_that("the bundled binary is used on Linux", {
    skip_if_not(Sys.info()[["sysname"]] == "Linux")
    withr::local_options(misha.bigWigToWig = NULL)
    expect_true(file.exists(get_bigWigToWig_bin("Linux")))
})

test_that("a converter on PATH is used off-Linux", {
    withr::local_options(misha.bigWigToWig = NULL)
    dir <- withr::local_tempdir()
    fake <- file.path(dir, "bigWigToWig")
    writeLines("#!/bin/sh\nexit 0", fake)
    Sys.chmod(fake, "0755")
    withr::local_envvar(PATH = dir)

    expect_equal(normalizePath(get_bigWigToWig_bin("Darwin")), normalizePath(fake))
})

test_that("off-Linux with no converter errors with install instructions", {
    withr::local_options(misha.bigWigToWig = NULL)
    withr::local_envvar(PATH = withr::local_tempdir()) # empty -> Sys.which finds nothing

    expect_error(get_bigWigToWig_bin("Darwin"), "cannot run on Darwin")
    expect_error(get_bigWigToWig_bin("Darwin"), "ucsc-bigwigtowig")
    expect_error(get_bigWigToWig_bin("Darwin"), "misha\\.bigWigToWig")
})
