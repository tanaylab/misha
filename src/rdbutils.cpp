/*
 * rdbutils.cpp
 *
 *  Created on: Mar 26, 2010
 *      Author: hoichman
 */

#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <new>
#include <setjmp.h>
#include <stdarg.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>

#ifndef _POSIX_C_SOURCE
	#define _POSIX_C_SOURCE 199309
	#include <time.h>
	#undef _POSIX_C_SOURCE
#endif

#include "rdbinterval.h"
#include "rdbutils.h"
#include "System.h"
#include "Thread.h"
#include "TGLException.h"

//#define CSTACK_DEFNS

#include <Rinterface.h>
#include <Rembedded.h>
#include <R_ext/Parse.h>

using namespace std;
using namespace rdb;
using namespace __gnu_cxx;

// A hack into R to change the default error report mechanism.
//
// Surely this an ugly hack! R_GlobalContext is an internal object managed by R. Its type is not even exposed. RCNTXT structure is defined
// in an internal header file (src/include/Defn.h) that we cannot access. We copied the definition of it from there and changed it a bit
// to compile.

typedef struct RCNTXT {
    struct RCNTXT *nextcontext;	/* The next context up the chain */
    int callflag;		/* The context "type" */
    sigjmp_buf cjmpbuf;		/* C stack and register information */
    int cstacktop;		/* Top of the pointer protection stack */
    int evaldepth;	        /* evaluation depth at inception */
    SEXP promargs;		/* Promises supplied to closure */
    SEXP callfun;		/* The closure called */
    SEXP sysparent;		/* environment the closure was called from */
    SEXP call;			/* The call that effected this context*/
    SEXP cloenv;		/* The environment */
    SEXP conexit;		/* Interpreted "on.exit" code */
    void (*cend)(void *);	/* C "on.exit" thunk */
    void *cenddata;		/* data for C "on.exit" thunk */
    void *vmax;		        /* top of R_alloc stack */
    int intsusp;                /* interrupts are suspended */
    SEXP handlerstack;          /* condition handler stack */
    SEXP restartstack;          /* stack of available restarts */
    struct RPRSTACK *prstack;   /* stack of pending promises */
#ifdef BYTECODE
    SEXP *nodestack;
# ifdef BC_INT_STACK
    IStackval *intstack;
# endif
#endif
    SEXP srcref;	        /* The source line in effect */
} RCNTXT, *context;

static const char *CHROM_FILE_PREFIX = "chr";
static const unsigned CHROM_FILE_PREFIX_LEN = strlen(CHROM_FILE_PREFIX);
const string rdb::TRACK_FILE_EXT = ".track";
const string rdb::INTERV_FILE_EXT = ".interv";

const int64_t        RdbInitializer::LAUNCH_DELAY = 50; // in msec
const int64_t        RdbInitializer::MEM_SYNC_DELAY = 100;
const int64_t        RdbInitializer::REPORT_INTERVAL_DELAY = 3000;
size_t               RdbInitializer::s_shm_size;
size_t               RdbInitializer::s_res_const_size;
size_t               RdbInitializer::s_res_var_size;
size_t               RdbInitializer::s_max_res_size;
size_t               RdbInitializer::s_max_mem_usage;
bool                 RdbInitializer::s_is_kid = false;
pid_t                RdbInitializer::s_parent_pid = 0;
sem_t                RdbInitializer::s_signal_sem;
sem_t               *RdbInitializer::s_shm_sem = SEM_FAILED;
sem_t               *RdbInitializer::s_alloc_suspend_sem = SEM_FAILED;
int                  RdbInitializer::s_kid_index;
vector<RdbInitializer::LiveStat>     RdbInitializer::s_running_pids;
vector<RdbInitializer::DeathStat>    RdbInitializer::s_dead_pids;
RdbInitializer::Shm *RdbInitializer::s_shm = (RdbInitializer::Shm *)MAP_FAILED;
struct sigaction     RdbInitializer::s_old_sigint_act;
struct sigaction     RdbInitializer::s_old_sigchld_act;
int                  RdbInitializer::s_ref_count = 0;
int                  RdbInitializer::s_sigint_fired = 0;
unsigned             RdbInitializer::s_protect_counter = 0;

// To delete all unreleased shared memory run command:
// ipcrm `ipcs -m | awk '$6==0 {print("-m", $2)}'`
RdbInitializer::RdbInitializer()
{
// disable R check stack limit: required if eval is called not from the main thread
//R_CStackLimit=-1;
	s_sigint_fired = 0;

	if (!s_ref_count) {
		m_old_umask = umask(07);

		s_shm_size = 0;
		s_is_kid = false;
		s_parent_pid = getpid();
		sem_init(&s_signal_sem, 0, 0);
		s_shm_sem = SEM_FAILED;
		s_alloc_suspend_sem = SEM_FAILED;
		s_shm = (Shm *)MAP_FAILED;
		s_kid_index = 0;
		s_running_pids.clear();
		s_dead_pids.clear();

		// These containers are altered in a signal handler of SIGCHLD.
		// Since malloc should be avoided in a signal handler let's reserve the memory in advance.
		s_running_pids.reserve(MAX_KIDS);
		s_dead_pids.reserve(MAX_KIDS);

		m_old_error_handler = TGLException::set_error_handler(TGLException::throw_error_handler);

		// install out of memory handler
		m_old_new_handler = set_new_handler(out_of_memory);

		struct sigaction new_act;

		// install a new SIGINT handler
		new_act.sa_handler = sigint_handler;
		sigemptyset(&new_act.sa_mask);
		new_act.sa_flags = SA_RESTART;
		sigaction(SIGINT, &new_act, &s_old_sigint_act);

		// install a new SIGCHLD handler
		new_act.sa_handler = sigchld_handler;
		sigemptyset(&new_act.sa_mask);
		new_act.sa_flags = SA_RESTART | SA_NOCLDSTOP;
		sigaction(SIGCHLD, &new_act, &s_old_sigchld_act);

		// record the currently opened file descriptors
		get_open_fds(m_old_open_fds);
	}

	// Default error message that error() function generates includes the caller ("Error in long-blalalalalalala: ...")
	// To prevent the caller to be printed we simply screw it up.
	// (This ugly hack was made after learning verrorcall_dflt() function in R source code (error.c).)
	RCNTXT *c = (RCNTXT *)R_GlobalContext;
	if (c) {
		c->call = R_NilValue;
		if (c->nextcontext)
			c->nextcontext->call = R_NilValue;
	}

	s_ref_count++;

	// deal with PROTECT / UNPROTECT
	m_old_protect_count = s_protect_counter;
}

RdbInitializer::~RdbInitializer()
{
	s_ref_count--;

	if (!s_ref_count) {
		// if this is a child, do not detach from shared memory and do not deallocate the semaphore:
		// if exception is thrown ~RdbInitializer is called first and then the child might need
		// to write the error into the shared memory
		if (!s_is_kid) {
			if (s_shm_sem != SEM_FAILED) {
				SemLocker sl(s_shm_sem);
				SigchldBlocker sb;

				// kill all the remaining child processes
				for (vector<LiveStat>::const_iterator ipid = s_running_pids.begin(); ipid != s_running_pids.end(); ++ipid) 
					kill(ipid->pid, SIGTERM);
			}

			// after SIGTERM is sent to all the kids let's wait till sigchld_hander() burries them all
			while (1) {
				{
					SigchldBlocker sb;
					if (s_running_pids.empty()) 
						break;
				}
				sem_wait(&s_signal_sem);
			}

			if (s_shm_sem != SEM_FAILED)
				sem_close(s_shm_sem); // semaphore should be already unlinked, only need to close it

			if (s_alloc_suspend_sem != SEM_FAILED) 
				sem_close(s_alloc_suspend_sem);

			if (s_shm != (Shm *)MAP_FAILED)
				munmap(s_shm, s_shm_size);
		}

		TGLException::set_error_handler(m_old_error_handler);

		// install old out of memory handler
		if (m_old_new_handler)
			set_new_handler(m_old_new_handler);

		// install old signal handlers
		sigaction(SIGINT, &s_old_sigint_act, NULL);
		sigaction(SIGCHLD, &s_old_sigchld_act, NULL);

		sem_destroy(&s_signal_sem);

		// close all file descriptors opened during the session
		set<int> open_fds;
		get_open_fds(open_fds);
		for (set<int>::const_iterator ifd = open_fds.begin(); ifd != open_fds.end(); ++ifd) {
			if (m_old_open_fds.find(*ifd) == m_old_open_fds.end())
				close(*ifd);
		}

		umask(m_old_umask);

		// do not revert to R's default error report
	}

	// deal with PROTECT / UNPROTECT
	unprotect(s_protect_counter - m_old_protect_count);
	s_protect_counter = m_old_protect_count;
}

string RdbInitializer::get_shm_sem_name()
{
	char buf[100];
	sprintf(buf, "misha-shm-%d", (int)getpid());
	return buf;
}

string RdbInitializer::get_alloc_suspend_sem_name()
{
	char buf[100];
	sprintf(buf, "misha-alloc-suspend-%d", (int)getpid());
	return buf;
}

void RdbInitializer::prepare4multitasking(size_t res_const_size, size_t res_var_size, size_t max_res_size, size_t max_mem_usage, unsigned num_planned_kids)
{
	if (num_planned_kids > MAX_KIDS) 
		verror("Too many child processes");

	if (s_shm_sem == SEM_FAILED) {
		sem_unlink(get_shm_sem_name().c_str()); // remove a semaphore if it was somehow not cleaned from the previous invocation of the lib
		if ((s_shm_sem = sem_open(get_shm_sem_name().c_str(), O_CREAT | O_EXCL, 0644, 1)) == SEM_FAILED)
			verror("sem_open failed: %s", strerror(errno));

		// Open a semaphore and right after that unlink it. The semaphore will be useful up until
		// the last process holding it (i.e. keeping it open) dies.
		sem_unlink(get_shm_sem_name().c_str());
	}

	if (s_alloc_suspend_sem == SEM_FAILED) {
		sem_unlink(get_alloc_suspend_sem_name().c_str()); // remove a semaphore if it was somehow not cleaned from the previous invocation of the lib
		if ((s_alloc_suspend_sem = sem_open(get_alloc_suspend_sem_name().c_str(), O_CREAT | O_EXCL, 0644, 0)) == SEM_FAILED)
			verror("sem_open failed: %s", strerror(errno));

		// Open a semaphore and right after that unlink it. The semaphore will be useful up until
		// the last process holding it (i.e. keeping it open) dies.
		sem_unlink(get_alloc_suspend_sem_name().c_str());
	}

	if (s_shm == (Shm *)MAP_FAILED) {
		s_res_const_size = res_const_size;
		s_res_var_size = res_var_size;
		s_max_res_size = max_res_size;
		s_max_mem_usage = max_mem_usage;
		s_shm_size = s_max_res_size + sizeof(Shm);

		s_shm = (Shm *)mmap(NULL, s_shm_size, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

		if (s_shm == (Shm *)MAP_FAILED) {
			if (errno == ENOMEM) 
				verror("Failed to allocate shared memory: %s\n"
					   "Memory usage of the library can be controlled via gmax.data.size option (see options, getOptions).",
					   strerror(errno));
			else
				verror("Failed to allocate shared memory: %s", strerror(errno));
		}

		s_shm->error_msg[0] = '\0';
		s_shm->res_offset = 0;
		s_shm->total_mem_usage = 0;
		s_shm->num_kids_running = num_planned_kids;
		s_shm->num_kids_suspended = 0;
		s_shm->untouchable_kid_idx = 0;

		for (int i = 0; i < MAX_KIDS; ++i) {
			s_shm->is_alive[i] = true;
			s_shm->mem_usage[i] = 0;
			s_shm->kid_progress[i] = 0;
			s_shm->kid_res_offset[i] = 0;
			s_shm->kid_res_num_records[i] = 0;
		}
	}
}

pid_t RdbInitializer::launch_process()
{
	if (s_shm_sem == SEM_FAILED || s_alloc_suspend_sem == SEM_FAILED || s_shm == (Shm *)MAP_FAILED)
		verror("Not prepared for multitasking");

	if (s_kid_index >= MAX_KIDS) 
		verror("Too many child processes");

	SigchldBlocker sb;

	check_interrupt();

	{
		SemLocker sl(s_shm_sem);
		if (s_shm->error_msg[0])
			verror("%s", s_shm->error_msg);
	}

	pid_t pid = fork(); 

	if (pid == -1)
		verror("fork failed: %s", strerror(errno));

	if (pid) { // a parent process
		s_running_pids.push_back(LiveStat(pid, s_kid_index));
		s_kid_index++;
	} else {   // a child process
		s_is_kid = true;

		sigaction(SIGINT, &s_old_sigint_act, NULL);
		sigaction(SIGCHLD, &s_old_sigchld_act, NULL);

		SEXP r_multitasking_stdout = GetOption(install("gmultitasking_stdout"), R_NilValue);

		if (!isLogical(r_multitasking_stdout) || !(int)LOGICAL(r_multitasking_stdout)[0])
			fclose(stdout);
		fclose(stderr);
		fclose(stdin);

		int64_t delta_mem_usage = get_unique_mem_usage(getpid()) - s_shm->mem_usage[s_kid_index];
		s_shm->mem_usage[s_kid_index] += delta_mem_usage;
		s_shm->total_mem_usage += delta_mem_usage;

		// set delay for child execution to prevent high initial bulk allocations
		struct timespec req;
		rdb::set_rel_timeout(LAUNCH_DELAY, req);
		for (int i = 0; i < s_kid_index; ++i) {
			nanosleep(&req, NULL);
			check_interrupt();
		}
	}

	return pid;
}

void RdbInitializer::wait_for_kids(rdb::IntervUtils &iu)
{
	int64_t delay_msec = LAUNCH_DELAY / 2;

	struct timespec timeout, last_progress_time, last_delay_change_time;
	int last_progress = -1;

	clock_gettime(CLOCK_REALTIME, &timeout);
	last_progress_time = last_delay_change_time = timeout;

	while (1) {
		clock_gettime(CLOCK_REALTIME, &timeout);
		set_abs_timeout(delay_msec, timeout);

		// sem_wait is interrupted by signal handlers
		while (sem_timedwait(&s_signal_sem, &timeout) < 0 && errno == EINTR)
			;

		check_interrupt();

		{
			SemLocker sl(s_shm_sem);
			if (s_shm->error_msg[0])
				verror("%s", s_shm->error_msg);
		}

		SigchldBlocker sb;
		while (s_dead_pids.size()) {
			const DeathStat &death_stat = s_dead_pids.back();

			{
				SemLocker sl(s_shm_sem);
				s_shm->num_kids_running--;
				if (death_stat.index >= 0) {
					s_shm->total_mem_usage -= s_shm->mem_usage[death_stat.index];
					s_shm->mem_usage[death_stat.index] = 0;
				}

			}
			if (WIFEXITED(death_stat.status))
				s_dead_pids.pop_back();
			else
				verror("Child process %d ended unexpectedly", (int)death_stat.pid);

			// increase the frequency of memory usage checks after child process death
			if (s_running_pids.size() > 1) {
				delay_msec = MEM_SYNC_DELAY / 2;
				clock_gettime(CLOCK_REALTIME, &last_delay_change_time);
			} else
				delay_msec = REPORT_INTERVAL_DELAY;
		}

		if (s_res_var_size) {
			// update of data_size is atomic => don't use a semaphore
			size_t res_num_records = 0;

			for (int i = 0; i < get_num_kids(); ++i)
				res_num_records += s_shm->kid_res_num_records[i];
			iu.verify_max_data_size(res_num_records, "Result");
		}

		if (s_running_pids.empty()) 
			break;

		if (s_running_pids.size() > 1) {
			int64_t total_mem_usage = s_shm->total_mem_usage;

			update_kids_mem_usage();

			// Check the memory consumption increase from the last measurement. If the memory limit is going to be
			// bridged till the next measurement => increase the monitoring rate.
			int64_t delta_mem_usage = s_shm->total_mem_usage - total_mem_usage;
			int64_t time2reach_limit = -1;

			if (delta_mem_usage > 0 && s_shm->total_mem_usage <= s_max_mem_usage)
				time2reach_limit = delay_msec * ((s_max_mem_usage - s_shm->total_mem_usage) / delta_mem_usage);

			if (time2reach_limit >= 0 && time2reach_limit < delay_msec) {
				delay_msec = min(MEM_SYNC_DELAY / 2, delay_msec);
				delay_msec = max(delay_msec, time2reach_limit);
				clock_gettime(CLOCK_REALTIME, &last_delay_change_time);
			} else if (delay_msec < REPORT_INTERVAL_DELAY && is_time_elapsed(2 * s_running_pids.size() * delay_msec, last_delay_change_time)) {
				// after (2 * s_running_pids.size() * initial_delay) increase the delay twice, i.e. lower down the frequency update_kids_mem_usage() is called
				// anyway do not drop the delay below what is required for progress reporting
				delay_msec = min(2 * delay_msec, REPORT_INTERVAL_DELAY);
				clock_gettime(CLOCK_REALTIME, &last_delay_change_time);
			}
		}

		if (is_time_elapsed(REPORT_INTERVAL_DELAY, last_progress_time)) {
			// update of progress is atomic => don't use a semaphore
			int progress = 0;

			for (int i = 0; i < get_num_kids(); ++i) 
				progress += s_shm->kid_progress[i];
			
			progress = progress / get_num_kids();
			if (progress < 100 && progress != last_progress)
				Rprintf("%d%%...", progress); 
			else {
				if (last_progress == -1) 
					Rprintf("0%%...");
				else
					Rprintf(".");
			}
			last_progress = progress;
			clock_gettime(CLOCK_REALTIME, &last_progress_time);
		}
		
		if (!s_shm->num_kids_running || s_shm->total_mem_usage < s_max_mem_usage) {
			// wake up suspended processes
			for (int i = 0; i < s_shm->num_kids_suspended; ++i) 
				sem_post(s_alloc_suspend_sem);
		}
	}

	if (last_progress >= 0) 
		Rprintf("100%%\n");
}

int64_t RdbInitializer::update_kids_mem_usage()
{
	// Update the memory usage of the kids with precise number.
	// Usually the kids maintain their own usage based on their own rough estimation.
	int64_t total_mem_usage = 0;

	for (vector<LiveStat>::const_iterator ipid = s_running_pids.begin(); ipid != s_running_pids.end(); ++ipid) {
		size_t mem_usage = get_unique_mem_usage(ipid->pid);

		if (mem_usage) {
			s_shm->mem_usage[ipid->index] = mem_usage;
//vdebug_print("\t\tUPDATE MEM USAGE OF %d to %ld\n", ipid->index, s_shm->mem_usage[ipid->index]);
			total_mem_usage += mem_usage;
		}
	}

	if (total_mem_usage) 
		s_shm->total_mem_usage = total_mem_usage;

	return total_mem_usage;
}

void RdbInitializer::handle_error(const char *msg)
{
	if (s_is_kid) {
		{
			SemLocker sl(s_shm_sem);
			if (!s_shm->error_msg[0]) { // write an error message only if there were no error messages before
				strncpy(s_shm->error_msg, msg, sizeof(s_shm->error_msg));
				s_shm->error_msg[sizeof(s_shm->error_msg)] = '\0';
			}
		}
		exit(1);
	} else
		error(msg);
}

void *RdbInitializer::allocate_res(size_t res_num_records)
{
	if (!s_is_kid)
		verror("allocate_res() cannot be called by parent process");

	update_res_data_size(res_num_records);

	SemLocker sl(s_shm_sem);
	s_shm->kid_res_offset[s_kid_index] = s_shm->res_offset;
	s_shm->res_offset += s_res_const_size + res_num_records * s_res_var_size;

	if (s_shm->res_offset > s_max_res_size) 
		verror("Result size exceeded the maximal allowed.\n"
				"Try to bound the scope of the function.\n"
				"Note: the maximum data size is controlled via gmax.data.size option (see options, getOptions).");

	return &s_shm->res + s_shm->kid_res_offset[s_kid_index];
}

void RdbInitializer::out_of_memory()
{
    if (!s_is_kid)
        Rprintf("\nOut of memory\n");
	verror("Out of memory\n");
}

void RdbInitializer::sigint_handler(int)
{
	++s_sigint_fired;

	// Normally this condition should be always true since the kid installs the default handler for SIGINT.
	// However due to race condition the old handler might still be in use.
	if (getpid() == s_parent_pid) {
		// wake up myself from wait_for_kids()
		sem_post(&s_signal_sem);
		printf("CTL-C!\n");
	}
}

void RdbInitializer::sigchld_handler(int)
{
	// Normally this condition should be always true since the kid installs the default handler for SIGINT.
	// However due to race condition the old handler might still be in use.
	if (getpid() == s_parent_pid) {
		int status;
		pid_t pid;

		while ((pid = waitpid((pid_t)-1, &status, WNOHANG)) > 0) {
			int kid_idx = -1;

			for (vector<LiveStat>::iterator ipid = s_running_pids.begin(); ipid != s_running_pids.end(); ++ipid) {
				if (ipid->pid == pid) {
					kid_idx = ipid->index;
					if (ipid != s_running_pids.end() - 1)
						*ipid = *(s_running_pids.end() - 1);
					s_running_pids.pop_back();

					// choose a new untouchable kid: the one with maximal memory consumption
					if (kid_idx == s_shm->untouchable_kid_idx && s_running_pids.size()) {
						int untouchable_kid_idx = s_running_pids.begin()->index;
						int64_t max_kid_mem_usage = s_shm->mem_usage[untouchable_kid_idx];
						for (vector<LiveStat>::iterator ipid2 = s_running_pids.begin() + 1; ipid2 < s_running_pids.end(); ++ipid2) {
							if (max_kid_mem_usage < s_shm->mem_usage[ipid2->index]) {
								untouchable_kid_idx = ipid->index;
								max_kid_mem_usage = s_shm->mem_usage[untouchable_kid_idx];
							}
						}
						s_shm->untouchable_kid_idx = untouchable_kid_idx;
					}
					break;
				}
			}

			if (kid_idx >= 0) {
				// malloc is forbidden in signal handlers, but the space was previously reserved for large number
				// of PIDs, therefore push_back does not initiate malloc
				s_dead_pids.push_back(DeathStat(pid, kid_idx, status));
				s_shm->is_alive[kid_idx] = false;
			}
		}

		sem_post(&s_signal_sem);
	}
}

void RdbInitializer::get_open_fds(set<int> &fds)
{
	DIR *dir = opendir("/proc/self/fd");
	struct dirent *dirp;

	fds.clear();
	while ((dirp = readdir(dir))) {
		char *endptr;
		int fd = strtol(dirp->d_name, &endptr, 10);
		if (!*endptr && fd != dirfd(dir)) // name is a number (it can be also ".", "..", whatever...)
			fds.insert(fd);
	}

	closedir(dir);
}

void RdbInitializer::vdebug_print(const char *fmt, ...)
{
	va_list ap;
	char buf[1000];

	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);

	SemLocker sl(s_shm_sem);
	printf("%s", buf);
}

void rdb::check_interrupt() throw (TGLException)
{
	if (RdbInitializer::s_sigint_fired)
		TGLError("Command interrupted!");
	monitor_memusage();
}

void rdb::rerror(const char *fmt, ...)
{
	va_list ap;
	char buf[1000];

	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);

	RdbInitializer::handle_error(buf);
}

void rdb::verror(const char *fmt, ...)
{
	va_list ap;
	char buf[1000];

	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);

	if (RdbInitializer::s_ref_count)
		TGLError("%s", buf);
	else
		RdbInitializer::handle_error(buf);
}

SEXP rdb::rprotect(SEXP &expr)
{
	if (expr != R_NilValue) {
		RdbInitializer::s_protect_counter++;
		return PROTECT(expr);
	}
	return expr;
}

void rdb::runprotect(int count)
{
	if (RdbInitializer::s_protect_counter < count)
		error("Number of calls to unprotect exceeds the number of calls to protect\n");
	UNPROTECT(count);
	RdbInitializer::s_protect_counter -= count;
}

void rdb::runprotect(SEXP &expr)
{
	if (expr != R_NilValue) {
		if (RdbInitializer::s_protect_counter < 1)
			error("Number of calls to unprotect exceeds the number of calls to protect\n");
		UNPROTECT_PTR(expr);
		expr = R_NilValue;
		RdbInitializer::s_protect_counter--;
	}
}

void rdb::runprotect(vector<SEXP> &exprs)
{
	for (vector<SEXP>::iterator iexpr = exprs.begin(); iexpr != exprs.end(); ++iexpr)
		runprotect(*iexpr);
}

void rdb::runprotect_all()
{
	if (RdbInitializer::s_protect_counter)
		UNPROTECT(RdbInitializer::s_protect_counter);
	RdbInitializer::s_protect_counter -= 0;
}

void rdb::get_chrom_files(const char *dirname, vector<string> &chrom_files)
{
	DIR *dir = opendir(dirname);

	if (!dir)
		verror("Failed to read directory %s: %s\n", dirname, strerror(errno));

	struct dirent *dirp;

	while ((dirp = readdir(dir))) {
		if (dirp->d_type == DT_REG && !strncmp(dirp->d_name, CHROM_FILE_PREFIX, CHROM_FILE_PREFIX_LEN))
			chrom_files.push_back(dirp->d_name);
	}

	closedir(dir);
}

const char *rdb::get_groot(SEXP envir)
{
	// no need to protect the returned value
	SEXP groot = findVar(install("GROOT"), envir);

	if (!isString(groot))
		verror("GROOT variable does not exist");

	return CHAR(STRING_ELT(groot, 0));
}

const char *rdb::get_gwd(SEXP envir)
{
	// no need to protect the returned value
	SEXP gwd = findVar(install("GWD"), envir);

	if (!isString(gwd))
		verror("GWD variable does not exist");

	return CHAR(STRING_ELT(gwd, 0));
}

const char *rdb::get_glib_dir(SEXP envir)
{
	// no need to protect the returned value
	SEXP glibdir = findVar(install(".GLIBDIR"), envir);

	if (!isString(glibdir))
		verror(".GLIBDIR variable does not exist");

	return CHAR(STRING_ELT(glibdir, 0));
}

string rdb::track2path(SEXP envir, const string &trackname)
{
	string path(trackname);
	for (string::iterator i = path.begin(); i != path.end(); ++i) {
		if (!is_R_var_char(*i))
			verror("Invalid track name %s. Only alphanumeric characters and _ are allowed in the name.", trackname.c_str());
		if (*i == '.')
			*i = '/';
	}
	return string(get_gwd(envir)) + "/" + path + TRACK_FILE_EXT;
}

string rdb::interv2path(SEXP envir, const string &intervname)
{
	string path(intervname);
	for (string::iterator i = path.begin(); i != path.end(); ++i) {
		if (!is_R_var_char(*i))
			verror("Invalid interval name %s. Only alphanumeric characters and _ are allowed in the name.", intervname.c_str());
		if (*i == '.')
			*i = '/';
	}
	return string(get_gwd(envir)) + "/" + path + INTERV_FILE_EXT;
}

string rdb::create_track_dir(SEXP envir, const string &trackname)
{
	string path = track2path(envir, trackname);
	if (mkdir(path.c_str(), 0777))
		verror("Cannot create track at %s: %s", path.c_str(), strerror(errno));

	return path;
}

string rdb::get_bounded_colname(const char *str, unsigned maxlen)
{
	string colname;

	maxlen = max(maxlen, 4u);
	if (strlen(str) > maxlen) {
		colname.assign(str, maxlen - 3);
		colname += "...";
	} else
		colname = str;
	return colname;
}

SEXP rdb::eval_in_R(SEXP parsed_command, SEXP envir)
{
	int check_error;
	SEXP res;

	rprotect(res = R_tryEval(parsed_command, envir, &check_error));
	if (check_error)
		verror(R_curErrorBuf());
	return res;
}

SEXP rdb::run_in_R(const char *command, SEXP envir)
{
	SEXP expr;
	SEXP parsed_expr;
	ParseStatus status;

	rprotect(expr = allocVector(STRSXP, 1));
	SET_STRING_ELT(expr, 0, mkChar(command));
	rprotect(parsed_expr = R_ParseVector(expr, -1, &status, R_NilValue));
	if (status != PARSE_OK)
		verror("Failed to parse expression \"%s\"", command);

	return eval_in_R(VECTOR_ELT(parsed_expr, 0), envir);
}

struct RSaneSerializeData {
	SEXP  rexp;
	FILE *fp;
};


static void RSaneSeserializeCallback(void *_data)
{
	RSaneSerializeData *data = (RSaneSerializeData *)_data;
	struct R_outpstream_st out;
	R_InitFileOutPStream(&out, data->fp, R_pstream_xdr_format, 2, NULL, NULL);
	R_Serialize(data->rexp, &out);
}

void rdb::RSaneSerialize(SEXP rexp, FILE *fp)
{
	RSaneSerializeData data;

	data.rexp = rexp;
	data.fp = fp;
    Rboolean ok = R_ToplevelExec(RSaneSeserializeCallback, &data);
	if (ok == FALSE)
		// We would like to print now the error contained in R_curErrorBuf(), however this error is automatically printed by R_ToplevelExec
		// and there's no way to prevent it without heavy hacking. On the other hand we want to abort the execution on error.
		// Solution: write a different error message. :)
		verror("Execution aborted");
}

void rdb::RSaneSerialize(SEXP rexp, const char *fname)
{
	FILE *fp = fopen(fname, "w");

	if (!fp)
		verror("Failed to open file %s: %s", fname, strerror(errno));

	RSaneSerialize(rexp, fp);
	fclose(fp);
}

struct RSaneUnserializeData {
	FILE *fp;
	SEXP  retv;
};

static void RSaneUnserializeCallback(void *_data)
{
	RSaneUnserializeData *data = (RSaneUnserializeData *)_data;
	struct R_inpstream_st in;
	R_InitFileInPStream(&in, data->fp, R_pstream_xdr_format, NULL, NULL);
	rprotect(data->retv = R_Unserialize(&in));
}

SEXP rdb::RSaneUnserialize(FILE *fp)
{
    RSaneUnserializeData data;

    data.fp = fp;
    data.retv = R_NilValue;

    Rboolean ok = R_ToplevelExec(RSaneUnserializeCallback, &data);
	if (ok == FALSE) 
		// We would like to print now the error contained in R_curErrorBuf(), however this error is automatically printed by R_ToplevelExec
		// and there's no way to prevent it without heavy hacking. On the other hand we want to abort the execution on error.
		// Solution: write a different error message. :)
		verror("Execution aborted");
	runprotect(1);
	return data.retv;
}

SEXP rdb::RSaneUnserialize(const char *fname)
{
	FILE *fp = fopen(fname, "r");

	if (!fp)
		verror("Failed to open file %s: %s", fname, strerror(errno));

	SEXP retv = RSaneUnserialize(fp);

	fclose(fp);
	return retv;
}

SEXP rdb::get_rvector_col(SEXP v, const char *colname, const char *varname, bool error_if_missing)
{
	SEXP colnames = getAttrib(v, R_NamesSymbol);

	if (!isVector(v) ||
		length(v) && (!isString(colnames) || length(colnames) != length(v)) ||
		!length(v) && !isNull(colnames))
		verror("Invalid format of %s", varname);

	int numcols = isNull(colnames) ? 0 : length(colnames);

	for (int i = 0; i < numcols; i++) {
		if (!strcmp(CHAR(STRING_ELT(colnames, i)), colname))
			return VECTOR_ELT(v, i);
	}

	if (error_if_missing)
		verror("Invalid format of %s: missing %s column", varname, colname);
	return R_NilValue;
}

