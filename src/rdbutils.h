
/*
 * rdbutils.h
 *
 *  Created on: Mar 26, 2010
 *      Author: hoichman
 */

#ifndef RDBUTILS_H_
#define RDBUTILS_H_

#if defined(__APPLE__)
#include <sys/time.h>
#include <mach/mach_time.h>
#include <mach/clock.h>
#include <mach/mach.h>

// Define clock identifiers for macOS if not already defined
#ifndef CLOCK_REALTIME
#define CLOCK_REALTIME 0
#endif

// Mark intentionally unused fallback helpers to suppress warnings in units
// that include this header but do not call them.
#if defined(__clang__) || defined(__GNUC__)
#define MISHA_MAYBE_UNUSED __attribute__((unused))
#else
#define MISHA_MAYBE_UNUSED
#endif

// Implement clock_gettime for macOS if needed (pre-Sierra 10.12)
#if !defined(HAVE_CLOCK_GETTIME)
static inline MISHA_MAYBE_UNUSED int clock_gettime(int clk_id, struct timespec *t)
{
	struct timeval tv;
	if (gettimeofday(&tv, NULL) < 0)
	{
		return -1;
	}
	t->tv_sec = tv.tv_sec;
	t->tv_nsec = tv.tv_usec * 1000;
	return 0;
}
#endif
#endif

#include <cstdint>
#include <string>
#include <vector>
#include <pthread.h>
#include <semaphore.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>

// Forward declare MultitaskingMode enum before including rdbinterval.h
// since rdbinterval.h needs to use this type
namespace rdb {
	// Multitasking execution modes
	enum MultitaskingMode : int {
		MT_MODE_SINGLE = 0,   // Single-threaded with unlimited heap
		MT_MODE_MMAP = 1      // Fast path with fixed mmap buffer (fork-based parallelism)
	};
}

#include "rdbinterval.h"
#include "Thread.h"

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>
#include <Rinterface.h>
#include <Rversion.h>

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

#define MISHA_PRENV(x) TAG(x)
#define MISHA_PRVALUE(x) CAR(x)
#define MISHA_PREXPR(x) R_BytecodeExpr(CDR(x))

#define LCONS1(a, b) Rf_lcons((a), (b))

#if R_VERSION < R_Version(4, 4, 1)
static inline SEXP Rf_allocLang(int n)
{
	if (n > 0)
		return LCONS1(R_NilValue, Rf_allocList(n - 1));
	else
		return R_NilValue;
}
#endif

#include "TGLException.h"

#define MISHA_EXIT_SIG SIGTERM

using namespace std;

//------------------------------- UTILITY FUNCTIONS -------------------------------------

namespace rdb {

class IntervUtils;

extern const string TRACK_FILE_EXT;
extern const string INTERV_FILE_EXT;

// should be used instead of R_CheckUserInterrupt. Throws exception if the command is interrupted.
void check_interrupt();

// adds timeout to the time that is already in req
void set_abs_timeout(int64_t delay_msec, struct timespec &req);

// sets timeout to req
void set_rel_timeout(int64_t delay_msec, struct timespec &req);

// returns true if current time exceeds start_time + delay
bool is_time_elapsed(int64_t delay_msec, const struct timespec &start_time);

// use rerror/verror instead of error!
//
// Neither returns, and neither can be made to return by anything a caller does. verror()
// inside a .Call throws (TGLError, which throws itself if the installed error handler
// hands control back); outside one it goes to handle_error(), which either ends a
// multitasking child through rexit() or reaches R through Rf_errorcall(). Every one of
// those is [[noreturn]] in turn, so the annotation here rests on construction rather than
// on the current contents of TGLException::s_error_handler - a mutable global with a
// public setter. Static analysers read it too: rchk otherwise walks the impossible
// fall-through past "runprotect(n); verror(...)" and reports the imbalance it invents.
[[noreturn]] void rerror(const char *fmt, ...);

[[noreturn]] void verror(const char *fmt, ...);

// Returns true the first time it is called with a given key within one top-level
// .Call, and always false inside a multitasking child process. Use it to run a
// diagnostic - and emit its warning - exactly once per user-level call: the track
// expression scanner, for instance, is built several times in the parent
// (prepare4multitasking, is_1d_iterator, the record estimate) and once more in
// every forked child.
bool once_per_call(const char *key);

// True inside a multitasking child process launched by distribute_task().
bool is_kid();

// Use rprotect instead of PROTECT!
SEXP rprotect(SEXP &expr);

// Like rprotect(SEXP&) but takes the object by value to avoid rchk
// "address taken" notes for local variables.
SEXP rprotect_ptr(SEXP expr);

// Unprotect the last "count" objects.
//
// SAFE ONLY WHILE NOTHING ELSE PROTECTS IN BETWEEN. The stack is shared with every
// other rprotect() in the call - IntervUtils' constructor pins ALLGENOME on it and
// never releases it, and convert_intervs() leaves its result on it - so "the n I
// protected" and "the last n" are the same objects only as long as nothing called in
// between protected anything of its own. When the number is conditional, or a helper
// runs in between, prefer runprotect(SEXP&) or protect_depth()/runprotect_to().
void runprotect(int count);

// Current depth of misha's protect stack, for use with runprotect_to().
unsigned protect_depth();

// Pops everything protected since protect_depth() returned "depth", and nothing below
// it. Unlike runprotect(int) the caller need not know how many objects the code in
// between chose to protect, so it stays correct when that changes. No-op if the stack
// is already at or below "depth".
void runprotect_to(unsigned depth);

// Unprotects object expr and sets it to R_NilValue. Works slower than runprotect(unsigned)!
void runprotect(SEXP &expr);

// Unprotects objects exprs and sets them to R_NilValue. Works slower than runprotect(unsigned)!
void runprotect(vector<SEXP> &exprs);

// Call runprotect_all if you wish to unprotect all object that are still protected
void runprotect_all();

struct SEXPCleaner {
    SEXPCleaner(SEXP &_var) : var(&_var) {}
    ~SEXPCleaner() { runprotect(*var); }
    SEXP *var;
};

void get_chrom_files(const char *dirname, vector<string> &chrom_files);

const char *get_groot(SEXP envir);

const char *get_gwd(SEXP envir);

const char *get_glib_dir(SEXP envir);

// Helper function to check if database is in indexed format
bool is_db_indexed(SEXP envir);

inline bool is_R_var_char(char c) { return isalnum(c) || c == '_' || c == '.'; }

// accepts track name, returns the path
string track2path(SEXP envir, const string &trackname);

// accepts track name, returns the path
string track2attrs_path(SEXP envir, const string &trackname);

// accepts track name, returns the path
string interv2path(SEXP envir, const string &intervname);

// Creates trackset and trackname directories using trackset.trackname. Verifies that the track name is valid.
string create_track_dir(SEXP envir, const string &trackname);

// the result is already protected
SEXP eval_in_R(SEXP parsed_command, SEXP envir);

// the result is already protected
SEXP run_in_R(const char *command, SEXP envir);

// This function writes R object to a file.
// Unlike R_Serialize function that just stops the execution if anything goes wrong (meaning: no clean up, destructors, etc.)
// RSaneSerialize throws an exception in case of error.
void RSaneSerialize(SEXP rexp, FILE *fp);
void RSaneSerialize(SEXP rexp, const char *fname);

// This function reads R object from a file. Object is expected to be saved using R's serialize() function or RSaneSerialize().
// THE RETURNED VALUE IS NOT PROTECTED. Protect it before the next allocation - all current
// callers either rprotect() it on the spot or store it into an already-protected object
// with nothing allocating in between. (The header used to claim the opposite; it could not
// be true, because R_ToplevelExec unwinds the protect stack to its own entry depth.)
// Unlike R_Unserialize function that just stops the execution if anything goes wrong (meaning: no clean up, destructors, etc.)
// RSaneUnserialize throws an exception in case of error.
SEXP RSaneUnserialize(FILE *fp);
SEXP RSaneUnserialize(const char *fname);

// ---------------------------------------------------------------------------
// The RSane* allocators, and when a raw R allocator is fine
// ---------------------------------------------------------------------------
//
// R's allocators raise an R error when they cannot satisfy a request, and an R
// error is a longjmp. A longjmp out of a .Call entry point skips every C++
// destructor between the throw point and R's top level - including
// ~RdbInitializer, which is what returns s_ref_count to its entry value and
// releases the shared memory, the child processes and the temporary files. A
// stuck s_ref_count is not a leak: this failure mode has already produced two
// defects in which a worker's exit through R's shutdown deleted the session
// tempdir, and the database inside it. This is invariant C3 of the lifetime
// contract: nothing inside an RdbInitializer scope may longjmp.
//
// Each RSane* wrapper runs the allocation inside R_ToplevelExec, so an
// allocation failure comes back as a TGLException instead. The entry point's
// own catch then unwinds the stack properly, running ~RdbInitializer on the
// way, and calls rerror() from a context where a longjmp is safe.
//
// WHEN TO USE ONE
//
//   Use an RSane* wrapper when the number of bytes the site requests grows
//   with the user's data - the length of a result column, one CHARSXP per
//   output row, a matrix sized by 4^k. Those are the requests that actually
//   run out of memory, and gseq.extract on a large interval set is the
//   canonical way a misha user does it.
//
// WHEN A RAW CALL IS FINE - and most of them are
//
//   1. Fixed-size requests. A 7-element column-name vector, Rf_mkString of a
//      literal, Rf_ScalarReal: R cannot be unable to find room for 8 bytes and
//      still be a working session. Wrapping these buys nothing and costs a
//      measured ~23 ns of R_ToplevelExec each.
//   2. Requests bounded by the shape of the database rather than by the call -
//      one element per chromosome, per column, per track. The widest genome in
//      the lab has 13,397 contigs; a STRSXP of chromosome names is under a
//      megabyte however large the query was.
//   3. Rf_mkChar(CHAR(STRING_ELT(src, i))) round-trips. R's CHARSXP cache
//      returns the *same* CHARSXP for a string that is already interned, and
//      every CHARSXP inside a live STRSXP is. Measured: 2,000,000 such calls
//      over a column of 26,000 distinct values allocated nothing at all.
//   4. Anywhere no RdbInitializer is live - C_revcomp, C_rev, C_comp,
//      C_intervals_coord_strings are pure R-to-R helpers, and a longjmp out of
//      one destroys nothing.
//
// Do not sweep. As of 5.11.22, 318 raw allocating calls remain in src/,
// and every one of them falls under one of those four heads. The worst defect
// found in the 2026-08 audit - a runprotect() count inside a loop that
// corrupted a database - was introduced by a commit whose stated purpose was
// hardening PROTECT usage across many sites at once.
//
// NOTE ON PROTECTION: like Rf_allocVector, none of these return a protected
// value. rprotect() the result before the next allocation.

// Replaces Rf_allocVector, which can fail on memory allocation and then R makes a longjmp, skipping all the destructors
SEXP RSaneAllocVector(SEXPTYPE type, R_xlen_t len);

// Same as above for Rf_mkChar, which allocates a CHARSXP and can fail the same way
SEXP RSaneMkChar(const char *str);

// Runs one R allocation under R_ToplevelExec and turns a failure into a
// TGLException. `alloc` must be a plain R allocator call and nothing else: it
// runs inside an R context, so it must not throw a C++ exception (that would
// unwind past endcontext) and must not call back into misha.
//
// This is the escape hatch for the allocators that have no named wrapper.
// There is deliberately no RSaneMkString / RSaneScalarReal / RSaneScalarInteger:
// every call site of those in misha allocates a fixed, tiny object, so naming
// a wrapper for them would only invite the mechanical sweep this comment
// exists to discourage.
template <typename Alloc>
SEXP RSaneAlloc(Alloc alloc)
{
	struct Payload { Alloc *alloc; SEXP retv; } payload = { &alloc, R_NilValue };

	Rboolean ok = R_ToplevelExec([](void *p) {
			Payload *d = (Payload *)p;
			d->retv = (*d->alloc)();
		}, &payload);

	if (!ok)
		verror("Allocation failed");
	return payload.retv;
}

// Same as RSaneAllocVector for Rf_mkCharLenCE. gseq.extract makes one of these
// per interval, each holding a whole DNA sequence - the largest allocations
// misha ever asks R for.
inline SEXP RSaneMkCharLenCE(const char *str, int len, cetype_t enc)
{
	return RSaneAlloc([&]() { return Rf_mkCharLenCE(str, len, enc); });
}

// Same as RSaneAllocVector for Rf_allocMatrix.
inline SEXP RSaneAllocMatrix(SEXPTYPE type, int nrow, int ncol)
{
	return RSaneAlloc([&]() { return Rf_allocMatrix(type, nrow, ncol); });
}

SEXP get_rvector_col(SEXP v, const char *colname, const char *varname, bool error_if_missing);

// Backward-compatible shims for R < 4.5.0
#if R_VERSION < R_Version(4, 5, 0)
static inline SEXP R_getVar(SEXP sym, SEXP rho, Rboolean inherits) {
    SEXP val = inherits ? Rf_findVar(sym, rho) : Rf_findVarInFrame(sym, rho);
    if (val == R_UnboundValue)
        Rf_error("object '%s' not found", CHAR(PRINTNAME(sym)));
    MARK_NOT_MUTABLE(val);
    return val;
}
static inline SEXP R_getVarEx(SEXP sym, SEXP rho, Rboolean inherits, SEXP ifnotfound) {
    SEXP val = inherits ? Rf_findVar(sym, rho) : Rf_findVarInFrame(sym, rho);
    if (val == R_UnboundValue)
        return ifnotfound;
    MARK_NOT_MUTABLE(val);
    return val;
}
#endif

// Helper: safely find a symbol in the package's .misha environment.
// Note: the returned SEXP is not protected. PROTECT it if you will perform
// any allocations before you are done using it.
static inline SEXP find_in_misha(SEXP envir, const char *name) {
    SEXP misha_env = R_NilValue;
    misha_env = rprotect_ptr(R_getVar(Rf_install(".misha"), envir, (Rboolean)TRUE));
    SEXP val = R_getVarEx(Rf_install(name), misha_env, (Rboolean)TRUE, R_UnboundValue);
    // rprotect_ptr() is a no-op on R_NilValue, so pop only what was actually pushed -
    // the same guard define_in_misha() below already carries. An unguarded
    // runprotect(1) would pop the caller's object instead.
    runprotect(misha_env != R_NilValue);
    return val;
}

// Helper: safely define a symbol in the package's .misha environment.
// Ensures both the target environment and the value are protected
// during the Rf_defineVar call.
static inline void define_in_misha(SEXP envir, const char *name, SEXP value) {
    SEXP misha_env = R_NilValue;
    misha_env = rprotect_ptr(R_getVar(Rf_install(".misha"), envir, (Rboolean)TRUE));
    SEXP tmp = value;
    tmp = rprotect_ptr(tmp);
    Rf_defineVar(Rf_install(name), tmp, misha_env);
    // rprotect_ptr() is a no-op on R_NilValue, so unprotect only what was actually
    // protected: defining a NULL value used to unprotect one frame too many, which
    // trips "Number of calls to unprotect exceeds the number of calls to protect".
    runprotect((misha_env != R_NilValue) + (tmp != R_NilValue));
}

// Helper: queue a diagnostic for .gcall() to raise once the .Call has returned.
//
// Neither Rf_warning() nor an R-level message() can be raised from here: a caller that
// catches either with an exiting handler (tryCatch, or options(warn = 2)) longjmps out
// of the C++ frame, skipping ~RdbInitializer and leaving misha's PROTECT counter
// inflated for the rest of the session.
//
// The queue is .misha$.GPENDING.DIAGNOSTICS: a list of length-2 character vectors,
// c(severity, text), in the order they were produced. `severity` is "message" or
// "warning". Appending rather than assigning one slot per severity means a second
// diagnostic within one call cannot silently overwrite the first, and a new severity
// costs nothing on either side of the boundary.
//
// A multitasking child gets nothing queued: the .misha environment it would write to is
// the fork's private copy, which dies with the child, so a diagnostic queued there can
// never be raised - the same reason once_per_call() returns false in a kid. Every site
// that can queue from a kid also runs in the parent (the parent converts the same
// intervals before it forks, and merges the kids' flags after), so the user still gets
// the diagnostic, once.
static inline void add_pending_diagnostic(SEXP envir, const char *severity, const char *text) {
    if (is_kid())
        return;

    SEXP prev = find_in_misha(envir, ".GPENDING.DIAGNOSTICS");
    int n = (prev != R_UnboundValue && TYPEOF(prev) == VECSXP) ? Rf_length(prev) : 0;

    // find_in_misha does not protect its result and the allocations below can collect
    SEXP held = R_NilValue;
    if (n)
        held = rprotect_ptr(prev);

    // RSaneAllocVector/RSaneMkChar rather than the plain R allocators: a failing
    // Rf_allocVector or Rf_mkChar longjmps out of here, which is precisely the unwind
    // this queue exists to avoid. These allocations are small, but the failure mode is
    // the same one, and it would strike while a diagnostic is being recorded.
    SEXP queue = R_NilValue;
    queue = rprotect_ptr(RSaneAllocVector(VECSXP, n + 1));
    for (int i = 0; i < n; ++i)
        SET_VECTOR_ELT(queue, i, VECTOR_ELT(held, i));

    SEXP entry = R_NilValue;
    entry = rprotect_ptr(RSaneAllocVector(STRSXP, 2));
    SET_STRING_ELT(entry, 0, RSaneMkChar(severity));
    SET_STRING_ELT(entry, 1, RSaneMkChar(text));
    SET_VECTOR_ELT(queue, n, entry);

    define_in_misha(envir, ".GPENDING.DIAGNOSTICS", queue);
    runprotect(n ? 3 : 2);
}

// Drops whatever add_pending_diagnostic() queued so far in this .Call. Used by the
// single-shard fallback (see rdb::SingleShard): the serial entry point re-runs the call
// from scratch and re-queues every diagnostic it warrants, so anything the multitasking
// entry point queued on its way to the fallback decision would otherwise be raised twice.
static inline void clear_pending_diagnostics(SEXP envir) {
    define_in_misha(envir, ".GPENDING.DIAGNOSTICS", R_NilValue);
}


void prepare4multitasking(uint64_t res_const_size, uint64_t res_var_size, uint64_t max_res_size, uint64_t max_mem_usage, unsigned num_planned_kids);

pid_t launch_process();

void wait_for_kids(IntervUtils &iu);

int get_num_kids();

void update_progress(unsigned char progress);

void update_res_data_size(uint64_t size);

// returns memory where the child process can write its result
void *allocate_res(uint64_t res_num_records);

// returns memory for the parent where the child process wrote its result
void *get_kid_res(int kid_index);

// returns result size of the child process in the number of data
uint64_t get_kid_res_size(int kid_index);

// keeps track of allocations in child processes; if the total memory consumption exceeds the limit,
// the child processes is suspended unless all the rest of the processes have been already suspended
void report_alloc(int64_t bytes);

// for child processes this function checks the memory usage and if needed suspends the child process up until
// the memory is freed
void monitor_memusage();

template<typename T> void pack_data(void *&ptr, const T &data, uint64_t n) {
	uint64_t size = sizeof(data) * n;
	memcpy(ptr, &data, size);
	ptr = (char *)ptr + size;
}

template<typename T> void unpack_data(void *&ptr, T &data, uint64_t n) {
	uint64_t size = sizeof(data) * n;
	memcpy(&data, ptr, size);
	ptr = (char *)ptr + size;
}

}

#define MAX_KIDS 1000
#define rreturn(retv) { if (RdbInitializer::is_kid()) rexit(); return(retv); }

namespace rdb {

// Thrown by a multitasking entry point that has decided not to multitask after all: either
// prepare4multitasking() planned a single shard, i.e. forking one kid would buy no
// parallelism, or the shared-memory arena could not be allocated and the retry ladder is
// exhausted. Catching it OUTSIDE the try block that owns the RdbInitializer is what makes
// the fallback safe: the throw unwinds that frame, so
// the RdbInitializer is destroyed before the serial entry point (with its own RdbInitializer)
// runs. Calling the serial entry point while the multitasking one still holds an
// RdbInitializer is not an option - an error inside it reaches R through Rf_error, whose
// longjmp skips the outer destructor and leaves s_ref_count stuck above zero. Every later
// multitasking call in that session then reuses stale shared memory and a stale kid index.
struct SingleShard {};

}

[[noreturn]] void rexit();

// Define RdbInitializer instance in your main function that is called by R.
// RdbInitializer should be defined inside "try-catch" statement that catches TGLException.
// RdbInitializer performs the following actions:
//   1. Installs a new SIGINT handler. ONE MUST CALL check_interrupt() INSTEAD OF R_CheckUserInterrupt()!!!!!!!
//   2. Installs out-of-memory handler.
//   3. suppresses the default error report behaviour.
//   4. Makes sure all file descriptors are closed on exit / error / interrupt.
//   5. Makes sure all objects are destructed on exit / error / interrupt.

class RdbInitializer {
public:
	RdbInitializer();
	~RdbInitializer();

	static bool   is_kid() { return s_is_kid; }
	static int    get_kid_idx() { return s_kid_index; }
    static void   get_open_fds(set<int> &fds);

	// The two counters that carry misha's lifetime invariant: both must be back at their
	// entry values when a .Call returns, and a longjmp past ~RdbInitializer is exactly what
	// leaves them stuck. Exposed so that C_lifetime_counters() (src/rdbutils.cpp) can report
	// them to R and a test can assert on the invariant itself rather than on a side effect
	// of it, such as the process umask.
	static int      get_ref_count() { return s_ref_count; }
	static unsigned get_protect_count() { return s_protect_counter; }

	// allows to safely write to stdout even from a child process
	// (before doing so please make sure launch_process() does not close stdout)
	static void vdebug_print(const char *fmt, ...);

private:
	struct LiveStat {
		pid_t pid;
		int   index;

		LiveStat(pid_t _pid, int _index) : pid(_pid), index(_index) {}
	};

	struct Shm {
		char          error_msg[10000];
		uint64_t        res_offset;
		int64_t       total_mem_usage;                 // cumulative memory usage of the kids
		uint64_t        num_kids_running;
		uint64_t        num_kids_suspended;
		int           untouchable_kid_idx;
		bool          is_alive[MAX_KIDS];
		int64_t       mem_usage[MAX_KIDS];
		unsigned char kid_progress[MAX_KIDS];          // progress report for each pid
		uint64_t        kid_res_offset[MAX_KIDS];        // offset for kid's result
		uint64_t        kid_res_num_records[MAX_KIDS];   // size of kid's result in number of data
		char          res;
	};

    struct SigBlocker {
        SigBlocker() {
            sigemptyset(&sigset);
            sigaddset(&sigset, SIGCHLD);
            sigaddset(&sigset, SIGINT);
            sigprocmask(SIG_BLOCK, &sigset, &oldsigset);
        }

        ~SigBlocker() { sigprocmask(SIG_UNBLOCK, &sigset, NULL); }

        sigset_t sigset;
        sigset_t oldsigset;
    };

	// all delays are in milliseconds
	static const int64_t        LAUNCH_DELAY;
	static const int64_t        MAX_LAUNCH_STAGGER;
	static const int64_t        MEM_SYNC_DELAY;
	static const int64_t        REPORT_INTERVAL_DELAY;

	static uint64_t               s_shm_size;
	static uint64_t               s_res_const_size;
	static uint64_t               s_res_var_size;
	static uint64_t               s_max_res_size;
	static uint64_t               s_max_mem_usage;
	static bool                 s_is_kid;
	static pid_t                s_parent_pid;
	static sem_t               *s_shm_sem;
	static sem_t               *s_alloc_suspend_sem;

	static int                  s_kid_index;
	static unsigned             s_num_planned_kids;
	static vector<LiveStat>     s_running_pids;
	static Shm                 *s_shm;

	static struct sigaction     s_old_sigint_act;
	static struct sigaction     s_old_sigchld_act;

	static int                  s_ref_count;
	static int                  s_sigint_fired;
	static unsigned             s_protect_counter;

	mode_t                      m_old_umask;
	TGLException::Error_handler m_old_error_handler;
	unsigned                    m_old_protect_count;
	set<int>                    m_old_open_fds;

	static string  get_shm_sem_name();
	static string  get_alloc_suspend_sem_name();
	static void    sigint_handler(int);
	static void    sigchld_handler(int);
	static void    prepare4multitasking(uint64_t res_const_size, uint64_t res_var_size, uint64_t max_res_size, uint64_t max_mem_usage, unsigned num_planned_kids);
	static pid_t   launch_process();
    static void    check_kids_state(bool ignore_errors);
	static void    wait_for_kids(rdb::IntervUtils &iu);
	static int64_t update_kids_mem_usage();
	static int     get_num_kids() { return s_kid_index; }
	[[noreturn]] static void handle_error(const char *msg);
	static void    update_progress(unsigned char progress);
	static void    update_res_data_size(uint64_t size);
	static void   *allocate_res(uint64_t res_num_records);
	static void   *get_kid_res(int kid_index);
	static uint64_t  get_kid_res_size(int kid_index);

	// report_alloc function keeps track of how much memory the child process has consumed so far.
	// Use positive value for new allocations and negative when the memory is freed.
	// If the total memory consumption of all the child processes exceeds the user defined limit, report_alloc
	// pauses the process (unless all the rest of the processes have been already paused).
	// This way the mechanism allows only one child to continue increasing the memory consumption.
	// When the child finishes and releases the memory, the paused processes are waken up.
	//
	// This mechanism should have been implemented using POSIX condition variables. However condition variables that can
	// be shared between processes are not compatible with some Linux systems. We don't want to rely on them. We implent the
	// mechanism with just one semaphore (s_alloc_suspend_sem) on which all the processes are going to sleep.
	// Since a check comes prior to the decision to sleep or allocate memory, without proper condition variables we are exposed
	// to race condition. But we consider the consequences of an error to be mild.
	// 1. If allocation is done whenever it should have been paused - not a big deal.
	// 2. If a child process is unnecesserily awaken - no problem: after a new check it will put itself to sleep again.
	// 3. If the process is paused due to race condition - this might create a deadlock. To battle it we are going to wake up
	//    everybody each 3 seconds in wait_for_kids(). So in the worst case even if we paused the process due to a very rare race condition,
	//    we will simply delay the execution by 3 seconds.
	//      Example: Consider 2 child processes and their interaction in time...
	//           T1 Child process 1 tries to make allocation and concludes that the limit is exceeded and child process 2 is still running.
	//           T2 Child process 2 dies
	//           T3 Parent process detects child process 2 death but doesn't wake up child process 1 because it is still not paused
	//           T4 Child process 1 pauses itself due to the decision it took at T1
	//           T5 (After 3 seconds) Parent process wakes up everybody (i.e. child process 1)
	//
	// Since report_alloc should practically be called before each and every memory allocation (or at least before big allocations)
	// incorporating this into the code might be somewhat problematic. We solve the problem by periodic checking of the memory usage
	// which is performed by the parent process. This check does not rely on what has been reported by report_alloc, but rather calls
	// Linux /proc interface to achieve the kids' memory consumption. Unfortunately this interface works only for Linux...
	// This check is performed every 3 seconds and the memory consumption counters are updated accordingly. By this one can guarantee
	// that high memory usage will be detected at some stage and the child processes will suspend themselves.
	//
	// The periodic memory consumption check has two issues:
	// 1. It works only on Linux. (We decided to ignore this problem for now.)
	// 2. It has a delay of 3 seconds. Since the memory consumption check is resource intensive we do not want to increase the rate of checking
	//    unless it is absolutely necessary.
	// 
	// The long delay between the checks creates an issue that child processes might breach the memory limit during these 3 seconds. The issue
	// is especially likely to happen right after the processes are created: very frequently a freshly created process performs significant
	// memory allocations that are required initiate its work. For example: each child process creates TrackExpressionScanner that might
	// load the whole track chromosome into memory (for example if the iterator is a sparse track). Thus if 10 processes are spawned,
	// 10 chromosomes might be loaded into memory which might cause high memory usage. Indeed after 3 seconds the parent will notice abnormal
	// memory usage, update memory usage counters and cause the child processes to pause themselves. Yet the damage is already done, the memory
	// is already consumed and the memory limit might be exceeded.
	//
	// We solve the problem via a "creation delay": before each child is created a delay is introduced. Child process number N delays itself
	// by N x creation_delay. Thus the first child is not being delayed at all while the last one delayes itself the longest.
	// By the end of the delay we hope that the previously created processes complete their initial bulk alloctions.
	// The memory consumption is then checked and the memory usage counters are updated accordingly. The newly created child in turn
	// will check the total memory usage and if it exceeds the limit it will pause itself. A similar delay is introduced after the child process
	// wakes up after suspension (which results from bridging the memory limit).
	//
	// The parent process in turn increases the rate of memory consumption checking while the processes are spawned or awaken. In the period
	// of "quiet" the rate slows down until it reaches 3 seconds.
	//
	// In addition the parent process chooses one process to be "untouchable", i.e. not suspendable. The untouchable child process never suspends
	// itself. This mechanism guarantees that each time only one child process is given a chance to finish its work without interrupts.
	// Switching from one process to another would not just possibly prolong the run-time, but it would certainly increase the total memory
	// consumption.
	//
	// After "untouchable" process dies, the next one is selected by choosing the process with the highest memory consumption.
	static void report_alloc(int64_t bytes);

	friend void rdb::check_interrupt();
	friend SEXP rdb::rprotect(SEXP &expr);
	friend void rdb::runprotect(int count);
	friend unsigned rdb::protect_depth();
	friend void rdb::runprotect(SEXP &expr);
	friend void rdb::runprotect(vector<SEXP> &exprs);
	friend void rdb::runprotect_all();
		friend SEXP rdb::rprotect_ptr(SEXP expr);
	friend void rdb::rerror(const char *fmt, ...);
	friend void rdb::verror(const char *fmt, ...);
	friend void rdb::prepare4multitasking(uint64_t res_const_size, uint64_t res_var_size, uint64_t max_res_size, uint64_t max_mem_usage, unsigned num_planned_kids);
	friend pid_t rdb::launch_process();
	friend void rdb::wait_for_kids(rdb::IntervUtils &iu);
	friend int rdb::get_num_kids();
	friend void rdb::update_progress(unsigned char progress);
	friend void rdb::update_res_data_size(uint64_t size);
	friend void *rdb::allocate_res(uint64_t res_num_records);
	friend void *rdb::get_kid_res(int kid_index);
	friend uint64_t rdb::get_kid_res_size(int kid_index);
	friend void rdb::report_alloc(int64_t bytes);

	friend class ChildShm;
};


// ------------------------------- IMPLEMENTATION --------------------------------

[[noreturn]] inline void rexit() {
	if (RdbInitializer::is_kid()){
		// Normally we should have called exit() here. However "R CMD check"
		// doesn't like calls to exit/abort/etc because they end R session
		// itself. It prints a Rf_warning message and packages with Rf_warning
		// messages cannot be submitted to CRAN. Yet the child process MUST end
		// the R sessions, that's the whole point. Solution? Send a signal to
		// itself. Fortunately "R CMD check" allows signals.
		//
		// IMPORTANT: Reset the signal handler to SIG_DFL before sending.
		// R installs its own SIGTERM handler that merely sets a flag instead
		// of terminating the process. Without this reset, the child survives
		// the signal and continues executing R-level code (e.g., gtrack.create
		// post-processing), causing file corruption and early parent return.
		struct sigaction sa;
		sa.sa_handler = SIG_DFL;
		sigemptyset(&sa.sa_mask);
		sa.sa_flags = 0;
		sigaction(MISHA_EXIT_SIG, &sa, NULL);

		// kill() returns 0 on success. POSIX only promises that the signal is delivered
		// before it returns if the signal is unblocked, so a blocked MISHA_EXIT_SIG - misha
		// never blocks it, but the mask a fork inherits is not misha's to assume - would
		// leave the child running R-level code on the parent's file descriptors, which is
		// the corruption the SIG_DFL reset above was added to stop. Unblock it first.
		sigset_t exit_sig;
		sigemptyset(&exit_sig);
		sigaddset(&exit_sig, MISHA_EXIT_SIG);
		sigprocmask(SIG_UNBLOCK, &exit_sig, NULL);
		kill(getpid(), MISHA_EXIT_SIG);

		// Unreachable in practice. It is here so that "this never returns" is a property of
		// the code rather than an argument about kill(): SIGKILL can be neither caught,
		// blocked nor ignored, and the loop leaves the compiler no path out of the branch.
		// A child that somehow cannot die then hangs - visible, and recoverable with Ctrl-C -
		// instead of falling through into R.
		for (;;)
			kill(getpid(), SIGKILL);
	}

	rdb::verror("rexit is called from parent process");
}

inline SEXP rdb::rprotect_ptr(SEXP expr)
{
    if (expr != R_NilValue) {
        RdbInitializer::s_protect_counter++;
        PROTECT(expr); // rchk: protect
        return expr;
    }
    return expr;
}

inline void rdb::set_abs_timeout(int64_t delay_msec, struct timespec &req)
{
	req.tv_nsec += delay_msec * 1000000L;
	req.tv_sec += req.tv_nsec / 1000000000L;
	req.tv_nsec %= 1000000000L;
}

inline void rdb::set_rel_timeout(int64_t delay_msec, struct timespec &req)
{
	req.tv_sec = delay_msec / 1000;
	req.tv_nsec = (delay_msec - req.tv_sec * 1000) * 1000000L;
}

inline bool rdb::is_time_elapsed(int64_t delay_msec, const struct timespec &start_time)
{
	struct timespec t1 = start_time;
	struct timespec t2;
	set_abs_timeout(delay_msec, t1);
	clock_gettime(CLOCK_REALTIME, &t2);
	return t2.tv_sec > t1.tv_sec || (t2.tv_sec == t1.tv_sec && t2.tv_nsec > t1.tv_nsec);
}

inline string rdb::track2attrs_path(SEXP envir, const string &trackname) {
	return rdb::track2path(envir, trackname) + "/.attributes";
}

inline void rdb::prepare4multitasking(uint64_t res_const_size, uint64_t res_var_size, uint64_t max_res_size, uint64_t max_mem_usage, unsigned num_planned_kids)
{
	RdbInitializer::prepare4multitasking(res_const_size, res_var_size, max_res_size, max_mem_usage, num_planned_kids);
}

inline pid_t rdb::launch_process() { return RdbInitializer::launch_process(); }

inline void rdb::wait_for_kids(IntervUtils &iu) { RdbInitializer::wait_for_kids(iu); }

inline int rdb::get_num_kids() { return RdbInitializer::get_num_kids(); }

inline void rdb::update_progress(unsigned char progress) { RdbInitializer::update_progress(progress); }

inline void rdb::update_res_data_size(uint64_t size) { RdbInitializer::update_res_data_size(size); }

inline void *rdb::allocate_res(uint64_t res_num_records) { return RdbInitializer::allocate_res(res_num_records); }

inline void *rdb::get_kid_res(int kid_index) { return RdbInitializer::get_kid_res(kid_index); }

inline uint64_t rdb::get_kid_res_size(int kid_index) { return RdbInitializer::get_kid_res_size(kid_index); }

inline void rdb::report_alloc(int64_t bytes) { RdbInitializer::report_alloc(bytes); }

inline void rdb::monitor_memusage() { report_alloc(0); }

inline void RdbInitializer::update_progress(unsigned char progress)
{
	if (s_is_kid)
		// update of progress is atomic => don't use a semaphore
		s_shm->kid_progress[s_kid_index] = progress;
}

inline void RdbInitializer::update_res_data_size(uint64_t size)
{
	if (s_is_kid)
		// update of progress is atomic => don't use a semaphore
		s_shm->kid_res_num_records[s_kid_index] = size;
}

inline void *RdbInitializer::get_kid_res(int kid_index)
{
	return &s_shm->res + s_shm->kid_res_offset[kid_index];
}

inline uint64_t RdbInitializer::get_kid_res_size(int kid_index)
{
	return s_shm->kid_res_num_records[kid_index];
}

inline void RdbInitializer::report_alloc(int64_t bytes)
{
	if (s_is_kid) {
//vdebug_print("%*s%d (%d): ATTEMPT TO ALLOC %ld, total: %ld, running: %ld, suspended: %ld\n", s_kid_index + 1, "", (int)s_kid_index, (int)getpid(), bytes,
//s_shm->total_mem_usage, s_shm->num_kids_running, s_shm->num_kids_suspended);
		if (s_kid_index != s_shm->untouchable_kid_idx) {  // never suspend untouchable kid
			while ((uint64_t)s_shm->total_mem_usage + bytes > s_max_mem_usage && s_shm->num_kids_running > 1) {
				{
					SemLocker sl(s_shm_sem);
					s_shm->num_kids_running--;
					s_shm->num_kids_suspended++;
				}
//vdebug_print("%*s%d (%d): SUSPENDING on ALLOC %ld, total: %ld, running: %ld, suspended: %ld\n", s_kid_index + 1, "", (int)s_kid_index, (int)getpid(), bytes,
//s_shm->total_mem_usage, s_shm->num_kids_running, s_shm->num_kids_suspended);

				while (sem_wait(s_alloc_suspend_sem) < 0 && errno == EINTR)
					;

				{
					SemLocker sl(s_shm_sem);
					s_shm->num_kids_suspended--;
					s_shm->num_kids_running++;
				}

				int num_preceding_kids = 0;
				for (int i = 0; i < s_kid_index; ++i) {
					if (s_shm->is_alive[i]) 
						++num_preceding_kids;
				}

				if (num_preceding_kids) {
					struct timespec req;
					rdb::set_rel_timeout(MEM_SYNC_DELAY, req);

					for (int i = 0; i < num_preceding_kids; ++i) {
						if (RdbInitializer::s_sigint_fired)
							TGLError("Command interrupted!");
						nanosleep(&req, NULL);
					}
				}

//vdebug_print("%*s%d (%d): WOKE up on ALLOC %ld, total: %ld, running: %ld, suspended: %ld\n", s_kid_index + 1, "", (int)s_kid_index, (int)getpid(), bytes,
//s_shm->total_mem_usage, s_shm->num_kids_running, s_shm->num_kids_suspended);
				if (RdbInitializer::s_sigint_fired)
					TGLError("Command interrupted!");
			}
		}

		// It would be too expensive to protect the next statements with a mutex/semaphore.
		// We're ready to suffer some errors in memory accounting on behalf of speed.
		if (bytes) { 
			s_shm->total_mem_usage += bytes;
			s_shm->mem_usage[s_kid_index] += bytes;
		}
//vdebug_print("%*s%d (%d): ALLOC %ld, total per process: %ld, total: %ld\n", s_kid_index + 1, "",
//(int)s_kid_index, (int)getpid(), bytes, s_shm->mem_usage[s_kid_index], s_shm->total_mem_usage);
	}
}

#endif /* RDBUTILS_H_ */

