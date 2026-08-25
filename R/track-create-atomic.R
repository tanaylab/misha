# Wrap a C++ track-create .gcall so the on-disk creation is atomic.
#
# Steps:
#  1. Compute a hidden tmp dir name: <parent>/.<basename>.tmp.<host>-<pid>.<rand>/
#     The leading "." and the embedded ".tmp." mean it is invisible to
#     gfind_tracks_n_intervals (which skips dirs with "." in the stem
#     before the .track suffix and treats *.track dirs as opaque).
#  2. Set .misha$.create_dir_override = tmp_dir so create_track_dir
#     mkdirs there instead of the final path. The C++ clears the slot
#     after reading it.
#  3. Run create_fn() (the .gcall to C++). C++ mkdirs tmp_dir and
#     writes per-chrom files into it.
#  4. On success: file.rename(tmp_dir, final_dir). Atomic within the
#     same parent directory on POSIX.
#  5. On failure: .gdb.trash(tmp_dir).
#
# Returns the final_dir path invisibly.
#
# Important: callers must do attribute writes, .gdb.add_track, and
# any post-create steps AFTER this function returns, so they operate
# on the final dir (which now exists after rename). Errors in those
# post-rename steps must call .gdb.trash(final_dir) themselves; this
# function only owns the create+rename atomicity.
#
# Concurrency: two sessions creating the same trackname will each
# write to a distinct tmp dir (host + PID + random suffix). The first to
# reach file.rename wins; the second sees the final_dir already
# present (rename onto a non-empty dir returns ENOTEMPTY/EEXIST on
# POSIX) and aborts with a "Refusing to overwrite existing" error,
# trashing its tmp.
#
# Orphan window: between file.rename and the caller's subsequent
# .gdb.add_track, the final dir exists on disk but is not in
# GTRACKS. Concurrent readers using the cached GTRACKS won't see
# it; a concurrent gdb.reload(rescan=TRUE) WILL pick it up. This
# is acceptable - the window is microseconds wide and no design
# without locking avoids it.
# Force `path` to stable storage. R has no fsync, so this is a thin .Call
# wrapper (src/GdbFsync.cpp). `recursive = TRUE` on a directory syncs every
# file under it and then the directory itself.
#
# Not swallowed: an fsync that fails means the writeback failed (EIO,
# ENOSPC on NFS), which is exactly the case where committing the rename
# would publish damaged data under the final name.
.gdb.fsync <- function(path, recursive = FALSE) {
    .gcall("gdb_fsync", path, isTRUE(recursive), .misha_env())
    invisible(NULL)
}

.gtrack.create_atomic <- function(trackname, create_fn) {
    final_dir <- .track_dir(trackname)
    parent <- dirname(final_dir)
    base <- basename(final_dir)
    if (!dir.exists(parent)) {
        stop(sprintf("Parent directory %s does not exist", parent), call. = FALSE)
    }
    # Clear leftovers from writers that died before they could clean up (see
    # .gdb.trash_sweep_old): a SIGKILL or a std::terminate never runs the
    # finally below, and the abandoned copy is invisible to gtrack.ls.
    try(.gdb.trash_sweep_old(parent), silent = TRUE)
    tmp_dir <- file.path(parent, .gdb.staging_name(base))

    assign(".create_dir_override", tmp_dir, envir = .misha)
    # Defensive: clear the slot on exit. Normally the C++ consumes it,
    # but if create_fn throws before any .gcall, the slot stays set.
    on.exit(
        {
            if (exists(".create_dir_override", envir = .misha, inherits = FALSE)) {
                assign(".create_dir_override", NULL, envir = .misha)
            }
        },
        add = TRUE
    )

    success <- FALSE
    tryCatch(
        {
            create_fn()
            if (!dir.exists(tmp_dir)) {
                stop(sprintf("Create completed but tmp dir %s missing", tmp_dir),
                    call. = FALSE
                )
            }
            if (dir.exists(final_dir)) {
                # .gconfirmtrackcreate should have caught this earlier;
                # treat as fatal.
                stop(sprintf("Refusing to overwrite existing %s", final_dir),
                    call. = FALSE
                )
            }
            # Durability, not just atomicity: rename() publishes the name
            # but guarantees nothing about the bytes reaching stable
            # storage, so a machine death (as opposed to a process death)
            # can leave the committed name pointing at a half-written
            # track. Sync the staged tree before the rename and the parent
            # directory after it.
            .gdb.fsync(tmp_dir, recursive = TRUE)
            if (!file.rename(tmp_dir, final_dir)) {
                stop(sprintf("Failed to rename %s -> %s", tmp_dir, final_dir),
                    call. = FALSE
                )
            }
            .gdb.fsync(parent)
            # Drop any stale cache entry for final_dir left over from a
            # previous track lifecycle at this path. Without this, readers
            # consulting GenomeTrack::s_index_cache may route to the wrong
            # layout (e.g., open a non-existent track.dat on what is now
            # a per-chrom track).
            .gdb.invalidate_dir_cache(final_dir)
            success <- TRUE
        },
        finally = {
            if (!success && dir.exists(tmp_dir)) {
                .gdb.trash(tmp_dir)
            }
        }
    )
    invisible(final_dir)
}
