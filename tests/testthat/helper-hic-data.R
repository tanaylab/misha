# Helper functions for creating and managing HiC test data
# These functions create synthetic 2D HiC tracks for testing based on
# real-world usage patterns from /home/nettam/projects/hic/lscripts/

#' Check if HiC test tracks exist in the test database
#'
#' @return logical vector indicating which tracks exist
check_hic_test_tracks <- function() {
    # Reload database to pick up any newly created tracks
    gdb.reload()

    required_tracks <- c(
        "hic.test_basic",
        "hic.test_score",
        "hic.test_shuffle",
        "insulation.test"
    )

    existing <- sapply(required_tracks, function(track) {
        tryCatch(
            {
                gtrack.exists(track)
            },
            error = function(e) FALSE
        )
    })

    names(existing) <- required_tracks
    existing
}

#' Check if test interval sets exist
#'
#' @return logical vector indicating which interval sets exist
check_hic_test_intervals <- function() {
    required_sets <- c("domains.test")

    existing <- sapply(required_sets, function(set) {
        tryCatch(
            {
                gintervals.exists(set)
            },
            error = function(e) FALSE
        )
    })

    names(existing) <- required_sets
    existing
}

#' Create synthetic HiC contact track with exponential distance decay
#'
#' This creates a 2D track with realistic distance-dependent contact frequency
#' modeled after real HiC data
#'
#' @param track_name Name of the track to create
#' @param chroms Vector of chromosome names (default: c("chr1", "chr2"))
#' @param chrom_sizes Named vector of chromosome sizes (default: 2Mb each)
#' @param resolution Bin resolution in bp (default: 10kb)
#' @param n_contacts_per_chr Number of contacts per chromosome (default: 5000)
#' @param decay_rate Exponential decay rate (default: -1.5)
#' @param min_dist Minimum distance between contacts (default: 1024)
create_hic_basic_track <- function(track_name = "hic.test_basic",
                                   chroms = c("chr1", "chr2"),
                                   chrom_sizes = c("chr1" = 2e6, "chr2" = 2e6),
                                   resolution = 1e4,
                                   n_contacts_per_chr = 5000,
                                   decay_rate = -1.5,
                                   min_dist = 1024) {
    message("Creating synthetic HiC track: ", track_name)

    # Generate contacts for each chromosome
    all_contacts <- list()

    for (chr in chroms) {
        chr_size <- chrom_sizes[chr]
        max_bins <- floor(chr_size / resolution)

        # Sample bin pairs with distance-dependent probability
        # Generate distances with exponential decay: more contacts at short distances
        max_dist <- chr_size
        # Use exponential with POSITIVE rate so we get more small values
        # rate = 1/mean_distance; higher rate = more short-distance contacts
        mean_dist <- 5e4 # Average contact distance ~50kb (creates strong short-range bias)
        distances <- rexp(n_contacts_per_chr * 10, rate = 1 / mean_dist)
        distances <- distances + min_dist # Shift to start at min_dist
        distances <- distances[distances >= min_dist & distances <= max_dist]
        distances <- distances[1:min(length(distances), n_contacts_per_chr)]

        # For each distance, randomly place the contact pair
        n_actual <- length(distances)
        starts1 <- sample(0:(chr_size - max(distances)), n_actual, replace = TRUE)
        starts2 <- starts1 + distances

        # Keep only valid contacts (within chromosome bounds)
        valid <- starts2 <= (chr_size - resolution)
        starts1 <- starts1[valid]
        starts2 <- starts2[valid]

        # Create data frame
        contacts <- data.frame(
            chrom1 = chr,
            start1 = floor(starts1 / resolution) * resolution,
            end1 = (floor(starts1 / resolution) + 1) * resolution,
            chrom2 = chr,
            start2 = floor(starts2 / resolution) * resolution,
            end2 = (floor(starts2 / resolution) + 1) * resolution
        )

        # Remove duplicates and ensure start1 <= start2
        contacts <- unique(contacts)
        swap <- contacts$start1 > contacts$start2
        if (any(swap)) {
            temp <- contacts$start1[swap]
            contacts$start1[swap] <- contacts$start2[swap]
            contacts$start2[swap] <- temp
            temp <- contacts$end1[swap]
            contacts$end1[swap] <- contacts$end2[swap]
            contacts$end2[swap] <- temp
        }

        all_contacts[[chr]] <- contacts
    }

    # Combine all contacts
    all_contacts_df <- do.call(rbind, all_contacts)
    row.names(all_contacts_df) <- NULL

    # Add contact values (count = 1 for each contact)
    values <- rep(1, nrow(all_contacts_df))

    # Create track from data frame
    gtrack.2d.create(track_name,
        description = "Synthetic HiC test track",
        intervals = all_contacts_df,
        values = values
    )

    message("Created ", nrow(all_contacts_df), " contacts in ", track_name)
    invisible(all_contacts_df)
}

#' Create synthetic score track from HiC track
#'
#' @param hic_track Source HiC track name
#' @param score_track Target score track name
#' @param mean_score Mean score value (default: 50)
#' @param sd_score Standard deviation of scores (default: 15)
create_hic_score_track <- function(hic_track = "hic.test_basic",
                                   score_track = "hic.test_score",
                                   mean_score = 50,
                                   sd_score = 15) {
    message("Creating score track: ", score_track)

    # Extract all intervals from HiC track
    intervals <- gextract(hic_track, intervals = gintervals.2d.all())

    if (nrow(intervals) == 0) {
        stop("No intervals found in ", hic_track)
    }

    # Generate scores
    intervals$score <- pmax(0, pmin(100, rnorm(nrow(intervals), mean_score, sd_score)))

    # Create track
    gtrack.2d.create(score_track,
        description = "Synthetic HiC score track",
        intervals = intervals[, c("chrom1", "start1", "end1", "chrom2", "start2", "end2")],
        values = intervals$score
    )

    message("Created score track with ", nrow(intervals), " scored intervals")
    invisible(intervals)
}

#' Create synthetic shuffled/expected track
#'
#' @param hic_track Source HiC track name
#' @param shuffle_track Target shuffled track name
#' @param noise_factor Multiplicative noise factor range (default: 0.8-1.2)
create_hic_shuffle_track <- function(hic_track = "hic.test_basic",
                                     shuffle_track = "hic.test_shuffle",
                                     noise_factor = c(0.8, 1.2)) {
    message("Creating shuffled track: ", shuffle_track)

    # Extract all intervals from HiC track
    intervals <- gextract(hic_track,
        intervals = gintervals.2d.all(),
        colnames = "contacts"
    )

    if (nrow(intervals) == 0) {
        stop("No intervals found in ", hic_track)
    }

    # Apply random noise
    intervals$shuffled <- intervals$contacts * runif(nrow(intervals),
        min = noise_factor[1],
        max = noise_factor[2]
    )

    # Create track
    gtrack.2d.create(shuffle_track,
        description = "Synthetic shuffled HiC track",
        intervals = intervals[, c("chrom1", "start1", "end1", "chrom2", "start2", "end2")],
        values = intervals$shuffled
    )

    message("Created shuffled track with ", nrow(intervals), " intervals")
    invisible(intervals)
}

#' Create synthetic 1D insulation track with periodic boundaries
#'
#' @param track_name Name of insulation track to create
#' @param chroms Chromosome names
#' @param chrom_sizes Named vector of chromosome sizes
#' @param resolution Bin resolution (default: 10kb)
#' @param n_domains_per_chr Number of domain boundaries per chromosome
#' @param boundary_strength Strength of boundaries (more negative = stronger)
create_insulation_track <- function(track_name = "insulation.test",
                                    chroms = c("chr1", "chr2"),
                                    chrom_sizes = c("chr1" = 2e6, "chr2" = 2e6),
                                    resolution = 1e4,
                                    n_domains_per_chr = 8,
                                    boundary_strength = -4) {
    message("Creating insulation track: ", track_name)

    all_bins <- list()

    for (chr in chroms) {
        chr_size <- chrom_sizes[chr]
        n_bins <- floor(chr_size / resolution)

        # Create bins
        starts <- seq(0, chr_size - resolution, by = resolution)
        ends <- starts + resolution

        # Initialize insulation values (baseline around -1)
        insulation <- rnorm(n_bins, mean = -1, sd = 0.3)

        # Add domain boundaries at regular intervals
        boundary_positions <- seq(1, n_bins, length.out = n_domains_per_chr + 1)
        boundary_positions <- unique(round(boundary_positions))

        # Set boundary values
        for (pos in boundary_positions) {
            if (pos > 0 && pos <= n_bins) {
                # Make boundaries stronger
                insulation[pos] <- boundary_strength + rnorm(1, 0, 0.2)
                # Smooth around boundaries
                if (pos > 1) insulation[pos - 1] <- (insulation[pos] - 1) / 2
                if (pos < n_bins) insulation[pos + 1] <- (insulation[pos] - 1) / 2
            }
        }

        bins <- data.frame(
            chrom = chr,
            start = starts[1:length(insulation)],
            end = ends[1:length(insulation)],
            insulation = insulation
        )

        all_bins[[chr]] <- bins
    }

    # Combine all bins
    all_bins_df <- do.call(rbind, all_bins)
    row.names(all_bins_df) <- NULL

    # Write to temporary file for import
    temp_file <- tempfile(fileext = ".txt")
    write.table(all_bins_df, temp_file,
        row.names = FALSE, col.names = FALSE,
        sep = "\t", quote = FALSE
    )

    # Import track
    gtrack.import(track_name,
        description = "Synthetic insulation track",
        file = temp_file,
        binsize = resolution
    )

    # Clean up temp file
    unlink(temp_file)

    message("Created insulation track with ", nrow(all_bins_df), " bins")
    invisible(all_bins_df)
}

#' Create synthetic domain intervals
#'
#' @param set_name Name of interval set to create
#' @param chroms Chromosome names
#' @param chrom_sizes Named vector of chromosome sizes
#' @param n_domains_per_chr Number of domains per chromosome
#' @param min_size Minimum domain size (default: 100kb)
#' @param max_size Maximum domain size (default: 500kb)
create_domains_intervals <- function(set_name = "domains.test",
                                     chroms = c("chr1", "chr2"),
                                     chrom_sizes = c("chr1" = 2e6, "chr2" = 2e6),
                                     n_domains_per_chr = 8,
                                     min_size = 1e5,
                                     max_size = 5e5) {
    message("Creating domain intervals: ", set_name)

    all_domains <- list()

    for (chr in chroms) {
        chr_size <- chrom_sizes[chr]

        # Generate domain sizes
        domain_sizes <- runif(n_domains_per_chr, min_size, max_size)

        # Place domains non-overlapping
        current_pos <- 0
        domains <- list()

        for (i in seq_along(domain_sizes)) {
            size <- domain_sizes[i]
            if (current_pos + size > chr_size) break

            domains[[i]] <- data.frame(
                chrom = chr,
                start = current_pos,
                end = current_pos + size
            )

            # Add gap between domains
            current_pos <- current_pos + size + runif(1, 0, min_size / 2)
        }

        if (length(domains) > 0) {
            all_domains[[chr]] <- do.call(rbind, domains)
        }
    }

    # Combine all domains
    all_domains_df <- do.call(rbind, all_domains)
    row.names(all_domains_df) <- NULL

    # Create directory for interval set if needed
    dir_name <- sub("\\..*", "", set_name) # Extract 'domains' from 'domains.test'
    tryCatch(gdir.create(dir_name, showWarnings = FALSE), error = function(e) NULL)

    # Save as interval set
    gintervals.save(set_name, all_domains_df)

    message("Created ", nrow(all_domains_df), " domains in ", set_name)
    invisible(all_domains_df)
}

#' Setup all HiC test tracks and interval sets
#'
#' This is the main function to call to ensure all test data exists
#'
#' @param force If TRUE, recreate tracks even if they exist
setup_hic_test_data <- function(force = FALSE) {
    message("Setting up HiC test data...")

    # Check which tracks exist
    existing_tracks <- check_hic_test_tracks()
    existing_sets <- check_hic_test_intervals()

    # Create tracks if needed
    if (force || !existing_tracks["hic.test_basic"]) {
        if (force && existing_tracks["hic.test_basic"]) {
            gtrack.rm("hic.test_basic", force = TRUE)
        }
        create_hic_basic_track()
    } else {
        message("hic.test_basic already exists, skipping")
    }

    if (force || !existing_tracks["hic.test_score"]) {
        create_hic_score_track()
    } else {
        message("hic.test_score already exists, skipping")
    }

    if (force || !existing_tracks["hic.test_shuffle"]) {
        create_hic_shuffle_track()
    } else {
        message("hic.test_shuffle already exists, skipping")
    }

    if (force || !existing_tracks["insulation.test"]) {
        create_insulation_track()
    } else {
        message("insulation.test already exists, skipping")
    }

    if (force || !existing_sets["domains.test"]) {
        create_domains_intervals()
    } else {
        message("domains.test already exists, skipping")
    }

    message("HiC test data setup complete!")

    # Return status
    list(
        tracks = check_hic_test_tracks(),
        intervals = check_hic_test_intervals()
    )
}

#' Cleanup all HiC test tracks and interval sets
#'
#' Use with caution - this will delete all test data
cleanup_hic_test_data <- function() {
    message("Cleaning up HiC test data...")

    tracks <- c("hic.test_basic", "hic.test_score", "hic.test_shuffle", "insulation.test")
    for (track in tracks) {
        if (gtrack.exists(track)) {
            gtrack.rm(track, force = TRUE)
            message("Removed track: ", track)
        }
    }

    sets <- c("domains.test")
    for (set in sets) {
        if (gintervals.exists(set)) {
            gintervals.rm(set, force = TRUE)
            message("Removed interval set: ", set)
        }
    }

    message("Cleanup complete!")
}

#' Get TAD domains from insulation track
#'
#' Extracts topologically associating domains (TADs) from an insulation track
#' by finding regions where insulation values are above a threshold or NA.
#'
#' @param insu_track Name of the insulation track
#' @param thresh Threshold value for insulation. Regions with values > thresh
#'   or NA are considered domains
#' @param iterator Bin size for the gscreen iterator (default: 500)
#'
#' @return A data frame of genomic intervals representing TAD domains
#'
#' @details
#' This function uses \code{gscreen} to identify regions where the insulation
#' score is either above the threshold or missing (NA). Higher insulation values
#' indicate weaker boundaries between TADs, so regions with high insulation
#' represent the interior of domains.
#'
#' @seealso \code{\link{gtrack.2d.get_insu_borders}}, \code{\link{gscreen}}
#'
#' @examples
#' \dontrun{
#' # Get domains from insulation track
#' domains <- gtrack.2d.get_insu_doms("insulation.test", thresh = -2)
#' }
#'
#' @noRd
gtrack.2d.get_insu_doms <- function(insu_track, thresh, iterator = 500) {
    doms <- gscreen(sprintf("is.na(%s) | %s > %f", insu_track, insu_track, thresh), iterator = iterator)
    return(doms)
}

#' Get TAD borders from insulation track
#'
#' Extracts topologically associating domain (TAD) borders from an insulation
#' track by finding regions where insulation values are below a threshold.
#'
#' @param insu_track Name of the insulation track
#' @param thresh Threshold value for insulation. Regions with values < thresh
#'   are considered domain borders
#' @param iterator Bin size for the gscreen iterator (default: 500)
#'
#' @return A data frame of genomic intervals representing TAD borders
#'
#' @details
#' This function uses \code{gscreen} to identify regions where the insulation
#' score is below the threshold and not missing. Lower (more negative) insulation
#' values indicate stronger boundaries between TADs, representing regions where
#' chromatin interactions are depleted.
#'
#' @seealso \code{\link{gtrack.2d.get_insu_doms}}, \code{\link{gscreen}}
#'
#' @examples
#' \dontrun{
#' # Get borders from insulation track
#' borders <- gtrack.2d.get_insu_borders("insulation.test", thresh = -2)
#' }
#'
#' @noRd
gtrack.2d.get_insu_borders <- function(insu_track, thresh, iterator = 500) {
    bords <- gscreen(sprintf("!is.na(%s) & %s < %f", insu_track, insu_track, thresh), iterator = iterator)
    return(bords)
}
