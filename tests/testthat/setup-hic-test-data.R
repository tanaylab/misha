# #!/usr/bin/env Rscript
# # Setup script to create test HiC tracks
# # Run this once to create the test data in the snapshot database

# # Load misha
# library(misha)

# # Set database root
# test_db <- "/net/mraid20/export/tgdata/db/tgdb/misha_snapshot/hg19"
# gsetroot(test_db)

# # Source helper functions (use correct path)
# helper_path <- if (file.exists("tests/testthat/helper-hic-data.R")) {
#     "tests/testthat/helper-hic-data.R"
# } else if (file.exists("helper-hic-data.R")) {
#     "helper-hic-data.R"
# } else {
#     stop("Cannot find helper-hic-data.R")
# }
# source(helper_path)

# # Create all test data
# message("Setting up HiC test tracks...")
# status <- tryCatch(
#     setup_hic_test_data(force = FALSE),
#     error = function(e) {
#         message("Error during setup: ", e$message)
#         list(tracks = check_hic_test_tracks(), intervals = check_hic_test_intervals())
#     }
# )

# # Print status
# message("\nTest data setup complete!")
# message("Tracks created:")
# print(status$tracks)
# message("\nInterval sets created:")
# print(status$intervals)

# message("\nYou can now run tests with: devtools::test(filter='2d-hic')")
