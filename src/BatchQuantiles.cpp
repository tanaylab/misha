// BatchQuantiles.cpp — batched multi-track genome-wide quantile scan.
//
// Uses BatchTrackScan<TopKQuantile> to scan N tracks in parallel threads.
// Phase 1 operates in fallback mode (full-vector storage + nth_element),
// matching the pre-refactor behavior bit-for-bit. Phase 2 adds top-K
// pruning + aggregator templating + intervals support.

#include "BatchTrackScan.h"
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
// TopKQuantile reducer. Phase 1: fallback-only (no top-K heap). State stores
// all non-NaN accepted values; finalize() runs nth_element per percentile.
// ---------------------------------------------------------------------------
struct TopKQuantile {
    struct Config {
        std::vector<double> percentiles;   // as provided (original order)
        uint32_t K = 0;
        bool use_fallback = true;          // Phase 1: always fallback
        bool top_side = true;
    };

    struct State {
        const Config *cfg = nullptr;
        uint64_t n_total = 0;
        std::vector<float> buf;            // full buffer in fallback mode

        void init(const Config &c, int /*chromid*/, int /*iterator_step*/) {
            cfg = &c;
            buf.reserve(1024);
        }

        void accept(float v, int64_t /*pos*/) {
            ++n_total;
            buf.push_back(v);
        }

        void boundary() {}

        bool prune(float /*upper*/, float /*lower*/) const { return false; }

        void merge(const State &o) {
            n_total += o.n_total;
            buf.insert(buf.end(), o.buf.begin(), o.buf.end());
        }
    };

    struct Result { std::vector<double> quantile_vals; };

    static constexpr bool needs_pruning = false;
    static constexpr bool needs_lower_bound = false;
};

static std::vector<double> topk_finalize(TopKQuantile::State &s)
{
    const auto &pctiles = s.cfg->percentiles;
    std::vector<double> out(pctiles.size(),
                            std::numeric_limits<double>::quiet_NaN());
    int64_t N = (int64_t)s.buf.size();
    if (N == 0) return out;
    for (size_t i = 0; i < pctiles.size(); ++i) {
        double p = pctiles[i];
        int64_t idx = (int64_t)std::floor(p * (double)(N - 1));
        if (idx < 0) idx = 0;
        if (idx >= N) idx = N - 1;
        std::nth_element(s.buf.begin(), s.buf.begin() + idx, s.buf.end());
        out[i] = (double)s.buf[idx];
    }
    return out;
}

extern "C" SEXP C_gquantiles_multi(
    SEXP _track_names,
    SEXP _percentiles,
    SEXP _iterator,
    SEXP _sshift,
    SEXP _eshift,
    SEXP _n_threads,
    SEXP _envir)
{
    try {
        RdbInitializer rdb_init;
        IntervUtils iu(_envir);

        int n_tracks = Rf_length(_track_names);
        std::vector<std::string> track_names(n_tracks);
        for (int i = 0; i < n_tracks; ++i)
            track_names[i] = CHAR(STRING_ELT(_track_names, i));

        int n_pctiles = Rf_length(_percentiles);
        std::vector<double> pctiles(REAL(_percentiles),
                                    REAL(_percentiles) + n_pctiles);

        int iterator_step = INTEGER(_iterator)[0];
        int sshift = INTEGER(_sshift)[0];
        int eshift = INTEGER(_eshift)[0];
        int n_threads_req = INTEGER(_n_threads)[0];

        if (iterator_step <= 0)
            verror("iterator must be a positive integer");

        const GenomeChromKey &chromkey = iu.get_chromkey();
        SEXP envir = iu.get_env();

        std::vector<std::string> track_dirs(n_tracks);
        std::vector<GenomeTrack::Type> track_types(n_tracks);
        for (int m = 0; m < n_tracks; ++m) {
            track_dirs[m] = track2path(envir, track_names[m]);
            track_types[m] = GenomeTrack::get_type(
                track_dirs[m].c_str(), chromkey, false);
        }

        std::vector<TopKQuantile::Config> configs(n_tracks);
        for (int m = 0; m < n_tracks; ++m) {
            configs[m].percentiles = pctiles;
            configs[m].K = 0;
            configs[m].use_fallback = true;
            configs[m].top_side = true;
        }

        unsigned hw = std::thread::hardware_concurrency();
        if (hw == 0) hw = 4;
        int n_threads;
        if (n_threads_req > 0) {
            n_threads = std::min(n_threads_req, n_tracks);
        } else {
            n_threads = (int)std::min((unsigned)n_tracks,
                                      std::min(hw, 40u));
        }
        if (n_threads < 1) n_threads = 1;

        ScanConfig scan;
        scan.func = WindowAggFunc::LSE;
        scan.iterator_step = iterator_step;
        scan.sshift = sshift;
        scan.eshift = eshift;
        scan.per_chrom_intervals = nullptr;
        scan.n_threads = n_threads;

        std::vector<BatchTrackScanTask<TopKQuantile>> tasks;
        run_batch_scan<TopKQuantile>(track_names, track_dirs, track_types,
                                     configs, scan, chromkey, tasks);

        for (auto &t : tasks) {
            if (!t.error_msg.empty())
                verror("Error processing track %s (chrom %d): %s",
                       t.track_name.c_str(), t.chromid, t.error_msg.c_str());
        }

        std::vector<TopKQuantile::State> per_track(n_tracks);
        for (int m = 0; m < n_tracks; ++m)
            per_track[m].init(configs[m], 0, iterator_step);
        for (auto &t : tasks)
            per_track[t.track_idx].merge(t.state);

        SEXP result;
        if (n_pctiles == 1)
            rprotect(result = Rf_allocVector(REALSXP, n_tracks));
        else
            rprotect(result = Rf_allocMatrix(REALSXP, n_tracks, n_pctiles));
        double *out = REAL(result);

        for (int m = 0; m < n_tracks; ++m) {
            auto qs = topk_finalize(per_track[m]);
            for (int p = 0; p < n_pctiles; ++p) {
                if (n_pctiles == 1)
                    out[m] = qs[p];
                else
                    out[m + (int64_t)p * n_tracks] = qs[p];
            }
        }

        if (n_pctiles == 1) {
            SEXP names;
            rprotect(names = Rf_allocVector(STRSXP, n_tracks));
            for (int i = 0; i < n_tracks; ++i)
                SET_STRING_ELT(names, i, STRING_ELT(_track_names, i));
            Rf_setAttrib(result, R_NamesSymbol, names);
        } else {
            SEXP rownames;
            rprotect(rownames = Rf_allocVector(STRSXP, n_tracks));
            for (int i = 0; i < n_tracks; ++i)
                SET_STRING_ELT(rownames, i, STRING_ELT(_track_names, i));

            SEXP colnames;
            rprotect(colnames = Rf_allocVector(STRSXP, n_pctiles));
            for (int i = 0; i < n_pctiles; ++i) {
                char buf[64];
                snprintf(buf, sizeof(buf), "%.6g", pctiles[i]);
                SET_STRING_ELT(colnames, i, Rf_mkChar(buf));
            }

            SEXP dimnames;
            rprotect(dimnames = Rf_allocVector(VECSXP, 2));
            SET_VECTOR_ELT(dimnames, 0, rownames);
            SET_VECTOR_ELT(dimnames, 1, colnames);
            Rf_setAttrib(result, R_DimNamesSymbol, dimnames);
        }

        return result;

    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const std::bad_alloc &) {
        rerror("Out of memory");
    } catch (const std::exception &e) {
        rerror("%s", e.what());
    }
    return R_NilValue;
}
