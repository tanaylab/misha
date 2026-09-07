#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "PottsParams.h"
#include "PottsModel.h"
#include "rdbutils.h"
#include "util.h" // log_sum_log

using namespace rdb;
using namespace std;

namespace {
enum PottsMode { PM_LSE = 0, PM_MAX = 1, PM_POS = 2, PM_COUNT = 3 };
}

extern "C" {

SEXP C_gseq_potts(SEXP r_seqs, SEXP r_params, SEXP r_mode, SEXP r_envir)
{
    try {
        RdbInitializer rdb_init;

        if (!Rf_isString(r_seqs))
            verror("gseq.potts: seqs must be a character vector");

        // The SAME parser the vtrack path uses, so the two cannot disagree
        // about what a model means.
        const PottsParams pp = PottsParams::parse(r_params, "gseq.potts");
        const PottsModel &model = pp.model;
        const PottsModel rc = model.rc();
        const int W = model.width();
        const int mode = Rf_asInteger(r_mode);
        const double thresh = pp.score_thresh;

        // bidirect wins over strand, as gseq.pwm does. When bidirect is FALSE
        // and strand is -1 only the reverse strand is read.
        const bool use_fwd = pp.bidirect || pp.strand_mode != -1;
        const bool use_rev = pp.bidirect || pp.strand_mode == -1;

        const R_xlen_t n = Rf_xlength(r_seqs);
        SEXP res = PROTECT(Rf_allocVector(REALSXP, n));
        double *out = REAL(res);

        vector<int8_t> codes;
        vector<int32_t> nbad;

        for (R_xlen_t r = 0; r < n; ++r) {
            SEXP el = STRING_ELT(r_seqs, r);
            if (el == NA_STRING) {
                out[r] = (mode == PM_COUNT) ? 0.0 : NA_REAL;
                continue;
            }
            const string target(CHAR(el));
            if ((int)target.size() < W) {
                out[r] = (mode == PM_COUNT) ? 0.0 : NA_REAL;
                continue;
            }
            potts_encode(target, codes, nbad);
            const size_t n_anchor = target.size() - (size_t)W + 1;

            double acc_lse = -numeric_limits<double>::infinity();
            bool have_lse = false;
            double best = -numeric_limits<double>::infinity();
            long long best_i = -1;
            int best_dir = 1;
            int count = 0;
            bool any = false;

            for (size_t i = 0; i < n_anchor; ++i) {
                // A Potts has no prior, so an ambiguous base leaves the anchor
                // unscorable rather than charged a fallback energy.
                if (nbad[i + (size_t)W] - nbad[i] != 0)
                    continue;
                const int8_t *c = &codes[i];
                double f = -numeric_limits<double>::infinity();
                double v = -numeric_limits<double>::infinity();
                if (use_fwd) f = model.score_codes(c);
                if (use_rev) v = rc.score_codes(c);

                // The strand union: log-sum-exp for lse/max/count, maximum for
                // pos, which has to name a strand. See the design doc's table -
                // inherited from the pwm family on purpose.
                double u;
                int dir = 1;
                if (use_fwd && use_rev) {
                    if (mode == PM_POS) {
                        u = (v > f) ? v : f;
                        dir = (v > f) ? -1 : 1;
                    } else {
                        u = f;
                        log_sum_log(u, v);
                    }
                } else if (use_fwd) {
                    u = f;
                } else {
                    u = v;
                    dir = -1;
                }

                any = true;
                if (mode == PM_LSE) {
                    // Seeded from the first scorable anchor rather than -inf:
                    // util.h's double log_sum_log() (util.h:57) has no isinf()
                    // guard, unlike the float overload at util.h:16, so
                    // log_sum_log(-inf, -inf) computes exp(NaN) = NaN. Every u
                    // here is finite by construction, but seeding this way
                    // means the accumulator can never touch that path.
                    if (!have_lse) {
                        acc_lse = u;
                        have_lse = true;
                    } else {
                        log_sum_log(acc_lse, u);
                    }
                }
                if (u > best) {
                    best = u;
                    best_i = (long long)i;
                    best_dir = dir;
                }
                if (mode == PM_COUNT && u >= thresh)
                    ++count;
            }

            if (mode == PM_COUNT)
                out[r] = (double)count;
            else if (!any)
                out[r] = NA_REAL;
            else if (mode == PM_LSE)
                out[r] = acc_lse;
            else if (mode == PM_MAX)
                out[r] = best;
            else // PM_POS: 1-based, signed by strand when bidirect
                out[r] = (double)(best_i + 1) * (pp.bidirect ? best_dir : 1);
        }

        UNPROTECT(1);
        return res;
    } catch (TGLException &e) {
        rerror("Error in C_gseq_potts: %s", e.msg());
    } catch (const std::exception &e) {
        rerror("Error in C_gseq_potts: %s", e.what());
    }
    rerror("Unknown error in C_gseq_potts");
    return R_NilValue;
}

} // extern "C"
