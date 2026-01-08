/*
 * TrackExpressionVars.h
 *
 *  Created on: Feb 9, 2012
 *      Author: hoichman
 */

#ifndef TRACKEXPRESSIONVARS_H_
#define TRACKEXPRESSIONVARS_H_

#include <cstdint>
#include <vector>
#include <string>
#include <string.h>
#include <memory>

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>

// Undefine R macros that conflict with C++ standard library
#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif
#ifdef warning
#undef warning
#endif

#include "BinFinder.h"
#include "Filter.h"
#include "FilterRegistry.h"
#include "GenomeTrack.h"
#include "GenomeTrackArrays.h"
#include "GenomeTrackInMemory.h"
#include "GInterval.h"
#include "GInterval2D.h"
#include "TrackExpressionIteratorBase.h"
#include "PWMScorer.h"
#include "KmerCounter.h"
#include "MaskedBpCounter.h"
#include "GenomeSeqFetch.h"
#include "rdbinterval.h"
#include "rdbutils.h"

// Forward declarations to avoid circular dependency
class TrackVarProcessor;
class IntervVarProcessor;
class ValueVarProcessor;
class SequenceVarProcessor;

using namespace std;

class TrackExpressionVars {
public:
    struct Binned_pv {
   		vector<double> bins;
        BinFinder  binfinder;
    };

    struct Iterator_modifier1D {
        enum Dimension { DIM_NONE, DIM1, DIM2 };

        Dimension dim;
        int64_t   sshift;
        int64_t   eshift;
        GInterval interval;
        bool      out_of_range;

        Iterator_modifier1D() : dim(DIM_NONE), sshift(0), eshift(0), out_of_range(false) {}
        bool operator==(const Iterator_modifier1D &o) const;

        void transform(const GInterval &interv, const GenomeChromKey &chromkey);
        void transform(const GInterval2D &interv, const GenomeChromKey &chromkey);
    };

    typedef vector<Iterator_modifier1D> Iterator_modifiers1D;

    struct Iterator_modifier2D {
        int64_t     sshift1;
        int64_t     eshift1;
        int64_t     sshift2;
        int64_t     eshift2;
        GInterval2D interval;
        bool        out_of_range;

        Iterator_modifier2D() : sshift1(0), eshift1(0), sshift2(0), eshift2(0), out_of_range(false) {}
        bool operator==(const Iterator_modifier2D &o) const;

        void transform(const GInterval2D &interv, const GenomeChromKey &chromkey);
    };

    typedef vector<Iterator_modifier2D> Iterator_modifiers2D;

    struct Track_n_imdf {
        string               name;
        GenomeTrack         *track{NULL};
        GenomeTrack::Type    type;
        vector<unsigned>     slice;
        GenomeTrackArrays::SliceFunctions slice_func;
        double               slice_percentile;
        Iterator_modifier1D *imdf1d{NULL};
        Iterator_modifier2D *imdf2d{NULL};
    };

    typedef vector<Track_n_imdf> Track_n_imdfs;

    struct Track_var {
        enum Val_func
        {
            REG,
            REG_MIN,
            REG_MAX,
            REG_NEAREST,
            STDDEV,
            SUM,
            QUANTILE,
            PV,
            PV_MIN,
            PV_MAX,
            WEIGHTED_SUM,
            OCCUPIED_AREA,
            PWM,
            PWM_MAX,
            PWM_MAX_POS,
            PWM_COUNT,
            KMER_COUNT,
            KMER_FRAC,
            MASKED_COUNT,
            MASKED_FRAC,
            MAX_POS_ABS,
            MAX_POS_REL,
            MIN_POS_ABS,
            MIN_POS_REL,
            EXISTS,
            SIZE,
            SAMPLE,
            SAMPLE_POS_ABS,
            SAMPLE_POS_REL,
            FIRST,
            FIRST_POS_ABS,
            FIRST_POS_REL,
            LAST,
            LAST_POS_ABS,
            LAST_POS_REL,
            NUM_FUNCS
        };

        static const char *FUNC_NAMES[NUM_FUNCS];

        string              var_name;
        SEXP                rvar{R_NilValue};
        double             *var;
        Val_func            val_func;
        double              percentile;
        bool                requires_pv;
        Binned_pv           pv_binned;
        Track_n_imdf       *track_n_imdf;
        std::unique_ptr<PWMScorer> pwm_scorer;
        std::unique_ptr<KmerCounter> kmer_counter;
        std::unique_ptr<MaskedBpCounter> masked_counter;
        char strand;
        // When the virtual track is sequence-based (PWM/KMER) we may still
        // have per-vtrack iterator modifiers; store them here.
        Iterator_modifier1D *seq_imdf1d{NULL};
        // Filter for masking genomic regions (applied after iterator modifiers)
        std::shared_ptr<Genome1DFilter> filter;
    };

    typedef vector<Track_var> Track_vars;

    struct Interv_var {
        enum Val_func { DIST, DIST_CENTER, DIST_EDGE, COVERAGE, NEIGHBOR_COUNT, NUM_FUNCS };

        static const char *FUNC_NAMES[NUM_FUNCS];

        string                     var_name;
        SEXP                       rvar{R_NilValue};
        double                    *var;
        Iterator_modifier1D       *imdf1d;
        Val_func                   val_func;

        GIntervals                 sintervs;
        GIntervals                 eintervs;
        GIntervals::const_iterator siinterv;
        GIntervals::const_iterator eiinterv;
        double                     dist_margin;
        // Filter for masking genomic regions (applied after iterator modifiers)
        std::shared_ptr<Genome1DFilter> filter;
    };

    typedef vector<Interv_var> Interv_vars;

    struct Value_var {
        enum Val_func {
            AVG, MIN, MAX, STDDEV, SUM, QUANTILE,
            NEAREST,
            EXISTS, SIZE,
            FIRST, LAST, SAMPLE,
            FIRST_POS_ABS, FIRST_POS_REL,
            LAST_POS_ABS, LAST_POS_REL,
            SAMPLE_POS_ABS, SAMPLE_POS_REL,
            MIN_POS_ABS, MIN_POS_REL,
            MAX_POS_ABS, MAX_POS_REL,
            NUM_FUNCS
        };

        static const char *FUNC_NAMES[NUM_FUNCS];

        string                     var_name;
        SEXP                       rvar{R_NilValue};
        double                    *var;
        Iterator_modifier1D       *imdf1d;
        Val_func                   val_func;
        double                     percentile;

        // Track object that holds the intervals and values
        std::shared_ptr<GenomeTrackInMemory> track;

        // Filter for masking genomic regions (applied after iterator modifiers)
        std::shared_ptr<Genome1DFilter> filter;
    };

    typedef vector<Value_var> Value_vars;

	TrackExpressionVars(rdb::IntervUtils &iu);
	~TrackExpressionVars();

	unsigned get_num_track_vars() const { return m_track_vars.size(); }
	unsigned get_num_interv_vars() const { return m_interv_vars.size(); }
	unsigned get_num_value_vars() const { return m_value_vars.size(); }

	const string &get_track_name(unsigned ivar) const { return m_track_vars[ivar].track_n_imdf->name; }
	GenomeTrack::Type get_track_type(unsigned ivar) const { return m_track_vars[ivar].track_n_imdf->type; }

	void parse_exprs(const vector<string> &track_exprs);
	void init(const TrackExpressionIteratorBase &expr_itr);
	void define_r_vars(unsigned size);
    const Track_var *var(const char *var_name) const;

    bool is_seq_variable(unsigned ivar) const;

    // Helper methods to check variable function types
    static bool is_sequence_based_function(Track_var::Val_func func);
    static bool is_pwm_function(Track_var::Val_func func);
    static bool is_kmer_function(Track_var::Val_func func);
    static bool is_masked_function(Track_var::Val_func func);

	void set_vars(const GInterval &interval, unsigned idx);
	void set_vars(const GInterval2D &interval, const DiagonalBand &band, unsigned idx);

private:
    rdb::IntervUtils &m_iu;
    string                  m_groot;
	Track_vars              m_track_vars;
	Interv_vars             m_interv_vars;
	Value_vars              m_value_vars;
	Track_n_imdfs           m_track_n_imdfs;
	Iterator_modifiers1D    m_imdfs1d;
	Iterator_modifiers2D    m_imdfs2d;
	GInterval               m_interval1d;
	GInterval2D             m_interval2d;
	DiagonalBand            m_band;

	// Shared sequence fetcher for all sequence-based vtracks to enable caching
	GenomeSeqFetch          m_shared_seqfetch;

	// Processors for different variable types (using pointers to avoid circular dependency)
	std::unique_ptr<TrackVarProcessor>       m_track_processor;
	std::unique_ptr<IntervVarProcessor>      m_interv_processor;
	std::unique_ptr<ValueVarProcessor>       m_value_processor;
	std::unique_ptr<SequenceVarProcessor>    m_sequence_processor;

	void                 parse_imdf(SEXP rvtrack, const string &vtrack, Iterator_modifier1D *imdf1d, Iterator_modifier2D *imdf2d);
	Iterator_modifier1D *add_imdf(const Iterator_modifier1D &imdf1d);
	Iterator_modifier2D *add_imdf(const Iterator_modifier2D &imdf2d);
	Track_n_imdf        &add_track_n_imdf(const string &track, GenomeTrack::Type track_type,
										  const vector<unsigned> &slice, GenomeTrackArrays::SliceFunctions slice_func, double slice_percentile,
										  const Iterator_modifier1D &imdf1d, const Iterator_modifier2D &imdf2d);
	Track_var           &add_track_var(const string &track);
	void                 add_vtrack_var(const string &track, SEXP rvtrack);
	Track_var           &add_vtrack_var_src_track(SEXP rvtrack, const string &vtrack, const string &track);
	Interv_var          &add_vtrack_var_src_interv(SEXP rvtrack, const string &vtrack, GIntervals &intervs1d, GIntervals2D &intervs2d);
	Value_var           &add_vtrack_var_src_value(SEXP rvtrack, const string &vtrack, GIntervals &intervs, vector<float> &vals);
	// Template function for attaching filter to any variable type with a 'filter' member
	template<typename VarType>
	void                 attach_filter_to_var(SEXP rvtrack, const string &vtrack, VarType &var);
	void                 register_track_functions();

	void start_chrom(const GInterval &interval);
	void start_chrom(const GInterval2D &interval);
	void set_vars(unsigned idx);


	bool is_var(const string &str, uint64_t start, uint64_t end) const { return (!start || !rdb::is_R_var_char(str[start - 1])) && (end == str.size() || !rdb::is_R_var_char(str[end])); }

    static int findListElementIndex(SEXP list, const char* name) {
        SEXP names = Rf_getAttrib(list, R_NamesSymbol);
        if (names == R_NilValue)
            rdb::verror("List must have named elements");
            
        int len = Rf_length(list);
        for (int i = 0; i < len; i++) {
            if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
                return i;
        }
        return -1;  // Element not found
    }
};


// -------------------------------------------------- IMPLEMENTATION -------------------------------------------------------

inline bool TrackExpressionVars::Iterator_modifier1D::operator==(const Iterator_modifier1D &o) const
{
	return dim == o.dim && sshift == o.sshift && eshift == o.eshift;
}

inline void TrackExpressionVars::Iterator_modifier1D::transform(const GInterval &interv, const GenomeChromKey &chromkey)
{
	interval.chromid = interv.chromid;
	interval.start = max(interv.start + sshift, (decltype(interv.start + sshift))0);
	interval.end = min(interv.end + eshift, (int64_t)chromkey.get_chrom_size(interv.chromid));
	interval.strand = interv.strand;
	out_of_range = interval.start >= interval.end;
}

inline void TrackExpressionVars::Iterator_modifier1D::transform(const GInterval2D &interv, const GenomeChromKey &chromkey)
{
	return dim == DIM1 ?
			transform(GInterval(interv.chromid1(), interv.start1(), interv.end1(), 0), chromkey) :
				transform(GInterval(interv.chromid2(), interv.start2(), interv.end2(), 0), chromkey);
}

inline bool TrackExpressionVars::Iterator_modifier2D::operator==(const Iterator_modifier2D &o) const
{
	return sshift1 == o.sshift1 && eshift1 == o.eshift1 && sshift2 == o.sshift2 && eshift2 == o.eshift2;
}

inline void TrackExpressionVars::Iterator_modifier2D::transform(const GInterval2D &interv, const GenomeChromKey &chromkey)
{
	int64_t start1 = max(interv.start1() + sshift1, (decltype(interv.start1() + sshift1))0);
	int64_t end1 = min(interv.end1() + eshift1, (int64_t)chromkey.get_chrom_size(interv.chromid1()));
	int64_t start2 = max(interv.start2() + sshift2, (decltype(interv.start2() + sshift2))0);
	int64_t end2 = min(interv.end2() + eshift2, (int64_t)chromkey.get_chrom_size(interv.chromid2()));
	interval.set(interv.chromid1(), start1, end1, interv.chromid2(), start2, end2);
	out_of_range = start1 >= end1 || start2 >= end2;
}

inline const TrackExpressionVars::Track_var *TrackExpressionVars::var(const char *var_name) const
{
    for (Track_vars::const_iterator ivar = m_track_vars.begin(); ivar != m_track_vars.end(); ivar++) {
        if (ivar->var_name == var_name)
            return &*ivar;
    }
    return NULL;
}

inline bool TrackExpressionVars::is_seq_variable(unsigned ivar) const {
    return m_track_vars[ivar].val_func == Track_var::PWM ||
           m_track_vars[ivar].val_func == Track_var::PWM_MAX ||
           m_track_vars[ivar].val_func == Track_var::PWM_MAX_POS ||
           m_track_vars[ivar].val_func == Track_var::PWM_COUNT ||
           m_track_vars[ivar].val_func == Track_var::KMER_COUNT ||
           m_track_vars[ivar].val_func == Track_var::KMER_FRAC ||
           m_track_vars[ivar].val_func == Track_var::MASKED_COUNT ||
           m_track_vars[ivar].val_func == Track_var::MASKED_FRAC;
}

// Helper methods to check variable function types
inline bool TrackExpressionVars::is_sequence_based_function(Track_var::Val_func func) {
    return func == Track_var::PWM || func == Track_var::PWM_MAX ||
           func == Track_var::PWM_MAX_POS || func == Track_var::PWM_COUNT ||
           func == Track_var::KMER_COUNT || func == Track_var::KMER_FRAC ||
           func == Track_var::MASKED_COUNT || func == Track_var::MASKED_FRAC;
}

inline bool TrackExpressionVars::is_pwm_function(Track_var::Val_func func) {
    return func == Track_var::PWM || func == Track_var::PWM_MAX ||
           func == Track_var::PWM_MAX_POS || func == Track_var::PWM_COUNT;
}

inline bool TrackExpressionVars::is_kmer_function(Track_var::Val_func func) {
    return func == Track_var::KMER_COUNT || func == Track_var::KMER_FRAC;
}

inline bool TrackExpressionVars::is_masked_function(Track_var::Val_func func) {
    return func == Track_var::MASKED_COUNT || func == Track_var::MASKED_FRAC;
}

#endif /* TRACKEXPRESSIONVARS_H_ */
