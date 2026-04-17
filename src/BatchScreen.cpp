// BatchScreen.cpp — batched multi-track threshold-screen scan.
//
// For each track, emits intervals where the windowed-aggregated value
// satisfies `value <op> threshold`. Consecutive passing scan positions
// within iterator-step distance merge into a single interval (same
// semantics as single-track gscreen). The main thread flattens all
// per-track intervals into a long data.frame with columns
// {chrom, start, end, track}.

#include "BatchTrackScan.h"
#include "rdbinterval.h"
#include "rdbutils.h"
#include "GenomeTrack.h"

#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <vector>

using namespace rdb;
using namespace batchscan;

// ---------------------------------------------------------------------------
// ThresholdScreen reducer. needs_pruning=true with bound-based skip for
// all operators except EQ; needs_lower_bound=true because LT/LE require
// window-min bounds.
// ---------------------------------------------------------------------------
struct ThresholdScreen {
    enum class CmpOp { LT = 0, LE = 1, EQ = 2, GE = 3, GT = 4 };

    struct Config { CmpOp op; float threshold; };

    struct State {
        const Config *cfg = nullptr;
        int chromid = -1;
        int iterator_step = 1;
        std::vector<GInterval> passing;   // flushed runs for this chrom
        int64_t cur_start = -1, cur_end = -1;

        void init(const Config &c, int cid, int it) {
            cfg = &c;
            chromid = cid;
            iterator_step = it;
            cur_start = cur_end = -1;
        }

        static bool compare(float v, float t, CmpOp op) {
            switch (op) {
                case CmpOp::LT: return v <  t;
                case CmpOp::LE: return v <= t;
                case CmpOp::EQ: return v == t;
                case CmpOp::GE: return v >= t;
                case CmpOp::GT: return v >  t;
            }
            return false;
        }

        void flush_cur() {
            if (cur_start >= 0) {
                // Each scan position represents a bin of width
                // iterator_step (matching the legacy gscreen bin semantics
                // where a passing position at pos emits an interval
                // [pos, pos + bin_size)). A run flushes as
                // [cur_start, cur_end + iterator_step) so consecutive
                // passing bins merge into contiguous intervals.
                passing.emplace_back(chromid, cur_start,
                                     cur_end + (int64_t)iterator_step,
                                     /*strand=*/(char)0);
                cur_start = cur_end = -1;
            }
        }

        void accept(float v, int64_t pos) {
            bool pass = compare(v, cfg->threshold, cfg->op);
            if (pass) {
                if (cur_start < 0 || pos > cur_end + iterator_step) {
                    flush_cur();
                    cur_start = pos;
                }
                cur_end = pos;
            } else {
                flush_cur();
            }
        }

        // Called by driver on interval-mask gap / chrom end. Critical for
        // correctness with fragmented intervals — without it, two passing
        // runs separated by a masked-out region would fuse into one
        // interval spanning positions that were never evaluated.
        void boundary() { flush_cur(); }

        // NaN aggregate → no value to compare against threshold; in legacy
        // gscreen semantics the comparison `NaN <op> threshold` is false
        // for every operator, which breaks any passing run.
        void nan_seen() { flush_cur(); }

        // EQ never prunes; GT/GE/LT/LE do. Returning false on EQ lets the
        // driver skip sliding-deque maintenance entirely for EQ predicates.
        bool pruning_active() const {
            return cfg && cfg->op != CmpOp::EQ;
        }

        bool prune(float upper, float lower) const {
            switch (cfg->op) {
                case CmpOp::GT: case CmpOp::GE: return upper <  cfg->threshold;
                case CmpOp::LT: case CmpOp::LE: return lower >  cfg->threshold;
                case CmpOp::EQ: return false;
            }
            return false;
        }

        // Pruned-but-valid positions must be counted in a reducer that
        // needs n_total for rank math; ThresholdScreen doesn't. No-op.
        void count_pruned() {}

        void merge(const State &o) {
            passing.insert(passing.end(), o.passing.begin(), o.passing.end());
        }
    };

    struct Result { std::vector<GInterval> intervals; };

    static constexpr bool needs_pruning = true;
    static constexpr bool needs_lower_bound = true;
};

static WindowAggFunc parse_func(SEXP _func)
{
    if (!Rf_isString(_func) || Rf_length(_func) != 1)
        verror("func must be a character scalar");
    const char *s = CHAR(STRING_ELT(_func, 0));
    if (!std::strcmp(s, "lse")) return WindowAggFunc::LSE;
    if (!std::strcmp(s, "avg")) return WindowAggFunc::AVG;
    if (!std::strcmp(s, "sum")) return WindowAggFunc::SUM;
    if (!std::strcmp(s, "max")) return WindowAggFunc::MAX;
    if (!std::strcmp(s, "min")) return WindowAggFunc::MIN;
    verror("func must be one of: lse, avg, sum, max, min");
    return WindowAggFunc::LSE;
}

// ---------------------------------------------------------------------------
// .Call entry: C_gscreen_multi
//
// Args (10):
//   1. _track_names    character vector (length N)
//   2. _iterator       integer scalar (bp step)
//   3. _sshift         integer scalar
//   4. _eshift         integer scalar
//   5. _n_threads      integer scalar (0 = auto)
//   6. _func           character scalar
//   7. _ops            integer vector of length N (CmpOp enum ints 0..4)
//   8. _thresholds     numeric vector of length N
//   9. _intervals      data.frame or NULL
//   10. _envir         R environment
//
// Returns: data.frame with columns chrom, start, end, track (long form).
// ---------------------------------------------------------------------------
extern "C" SEXP C_gscreen_multi(
    SEXP _track_names, SEXP _iterator, SEXP _sshift, SEXP _eshift,
    SEXP _n_threads, SEXP _func, SEXP _ops, SEXP _thresholds,
    SEXP _intervals, SEXP _envir)
{
    try {
        RdbInitializer rdb_init;
        IntervUtils iu(_envir);

        int n_tracks = Rf_length(_track_names);
        if (Rf_length(_ops) != n_tracks)
            verror("ops must have the same length as track_names");
        if (Rf_length(_thresholds) != n_tracks)
            verror("thresholds must have the same length as track_names");

        std::vector<std::string> track_names(n_tracks);
        for (int i = 0; i < n_tracks; ++i)
            track_names[i] = CHAR(STRING_ELT(_track_names, i));

        int iterator_step = INTEGER(_iterator)[0];
        int sshift = INTEGER(_sshift)[0];
        int eshift = INTEGER(_eshift)[0];
        int n_threads_req = INTEGER(_n_threads)[0];
        if (iterator_step <= 0)
            verror("iterator must be a positive integer");

        WindowAggFunc func = parse_func(_func);

        std::vector<ThresholdScreen::Config> configs(n_tracks);
        for (int m = 0; m < n_tracks; ++m) {
            int op_i = INTEGER(_ops)[m];
            if (op_i < 0 || op_i > 4)
                verror("ops[%d] = %d not in 0..4", m + 1, op_i);
            configs[m].op = (ThresholdScreen::CmpOp)op_i;
            configs[m].threshold = (float)REAL(_thresholds)[m];
        }

        const GenomeChromKey &chromkey = iu.get_chromkey();
        SEXP envir = iu.get_env();

        std::vector<std::string> track_dirs(n_tracks);
        std::vector<GenomeTrack::Type> track_types(n_tracks);
        for (int m = 0; m < n_tracks; ++m) {
            track_dirs[m] = track2path(envir, track_names[m]);
            track_types[m] = GenomeTrack::get_type(
                track_dirs[m].c_str(), chromkey, false);
        }

        const int n_chroms = (int)chromkey.get_num_chroms();
        std::vector<std::vector<GInterval>> per_chrom(n_chroms);
        bool use_intervals = !Rf_isNull(_intervals);
        if (use_intervals) {
            GIntervalsFetcher1D *i1d = nullptr;
            GIntervalsFetcher2D *i2d = nullptr;
            iu.convert_rintervs(_intervals, &i1d, &i2d);
            std::unique_ptr<GIntervalsFetcher1D> g1(i1d);
            std::unique_ptr<GIntervalsFetcher2D> g2(i2d);
            if (i2d && i2d->size() > 0)
                verror("intervals must be 1D");
            i1d->sort();
            i1d->unify_overlaps();
            for (int c = 0; c < n_chroms; ++c) {
                if (i1d->size(c) == 0) continue;
                i1d->begin_chrom_iter(c);
                for (auto it = i1d->get_chrom_begin();
                     it != i1d->get_chrom_end(); ++it) {
                    per_chrom[c].push_back(*it);
                }
            }
        }

        unsigned hw = std::thread::hardware_concurrency();
        if (hw == 0) hw = 4;
        int n_threads;
        if (n_threads_req > 0) {
            n_threads = std::min(n_threads_req, n_tracks);
        } else {
            n_threads = (int)std::min((unsigned)n_tracks, std::min(hw, 40u));
        }
        if (n_threads < 1) n_threads = 1;

        ScanConfig scan;
        scan.func = func;
        scan.iterator_step = iterator_step;
        scan.sshift = sshift;
        scan.eshift = eshift;
        scan.per_chrom_intervals = use_intervals ? &per_chrom : nullptr;
        scan.n_threads = n_threads;

        BatchTrackScanResult<ThresholdScreen> scan_result;
        run_batch_scan<ThresholdScreen>(track_names, track_dirs, track_types,
                                        configs, scan, chromkey, scan_result);

        for (int m = 0; m < n_tracks; ++m) {
            if (!scan_result.error_messages[m].empty())
                verror("Error processing track %s: %s",
                       track_names[m].c_str(),
                       scan_result.error_messages[m].c_str());
        }

        // Flatten per-track intervals into a long vector.
        int64_t total_intervals = 0;
        for (auto &st : scan_result.per_track_states)
            total_intervals += (int64_t)st.passing.size();

        // Return a track_idx column (0-based into the input vector) instead
        // of a track-name column. The R wrapper uses this to map each row
        // back to the caller's original expression string — unambiguous even
        // when multiple expressions share the same underlying source track.
        SEXP chrom_col, start_col, end_col, track_idx_col;
        rprotect(chrom_col = Rf_allocVector(STRSXP, total_intervals));
        rprotect(start_col = Rf_allocVector(REALSXP, total_intervals));
        rprotect(end_col = Rf_allocVector(REALSXP, total_intervals));
        rprotect(track_idx_col = Rf_allocVector(INTSXP, total_intervals));

        int64_t row = 0;
        for (int m = 0; m < n_tracks; ++m) {
            auto &ivs = scan_result.per_track_states[m].passing;
            for (auto &g : ivs) {
                SET_STRING_ELT(chrom_col, row,
                               Rf_mkChar(iu.id2chrom(g.chromid).c_str()));
                REAL(start_col)[row] = (double)g.start;
                REAL(end_col)[row] = (double)g.end;
                INTEGER(track_idx_col)[row] = m;
                ++row;
            }
        }

        // Build data.frame (list of 4 columns with class "data.frame" and
        // row.names attribute).
        SEXP df;
        rprotect(df = Rf_allocVector(VECSXP, 4));
        SET_VECTOR_ELT(df, 0, chrom_col);
        SET_VECTOR_ELT(df, 1, start_col);
        SET_VECTOR_ELT(df, 2, end_col);
        SET_VECTOR_ELT(df, 3, track_idx_col);

        SEXP names;
        rprotect(names = Rf_allocVector(STRSXP, 4));
        SET_STRING_ELT(names, 0, Rf_mkChar("chrom"));
        SET_STRING_ELT(names, 1, Rf_mkChar("start"));
        SET_STRING_ELT(names, 2, Rf_mkChar("end"));
        SET_STRING_ELT(names, 3, Rf_mkChar("track_idx"));
        Rf_setAttrib(df, R_NamesSymbol, names);

        SEXP rownames;
        rprotect(rownames = Rf_allocVector(INTSXP, 2));
        INTEGER(rownames)[0] = NA_INTEGER;
        INTEGER(rownames)[1] = -(int)total_intervals;
        Rf_setAttrib(df, R_RowNamesSymbol, rownames);

        SEXP cls;
        rprotect(cls = Rf_allocVector(STRSXP, 1));
        SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
        Rf_setAttrib(df, R_ClassSymbol, cls);

        return df;

    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const std::bad_alloc &) {
        rerror("Out of memory");
    } catch (const std::exception &e) {
        rerror("%s", e.what());
    }
    return R_NilValue;
}
