test_that(".gtf_to_genepred_path() honours MISHA_GTF_TO_GENEPRED env var", {
    stub <- tempfile(fileext = ".sh")
    writeLines("#!/bin/sh\necho hi", stub)
    Sys.chmod(stub, "0755")
    withr::with_envvar(c(MISHA_GTF_TO_GENEPRED = stub), {
        expect_equal(normalizePath(.gtf_to_genepred_path()), normalizePath(stub))
    })
})

test_that(".gtf_to_genepred_path() errors if env points at non-existent file", {
    withr::with_envvar(c(MISHA_GTF_TO_GENEPRED = "/no/such/path"), {
        expect_error(.gtf_to_genepred_path(), "non-existent")
    })
})

test_that(".gtf_to_genepred_path() returns NULL when no env and no cache", {
    skip_if_not_installed("withr")
    # Use a temporary HOME so R_user_dir resolves to a clean cache dir.
    tmp_home <- tempfile()
    dir.create(tmp_home)
    withr::with_envvar(
        c(
            MISHA_GTF_TO_GENEPRED = "",
            R_USER_CACHE_DIR = tmp_home,
            XDG_CACHE_HOME = tmp_home,
            HOME = tmp_home
        ),
        {
            expect_null(.gtf_to_genepred_path())
        }
    )
})
