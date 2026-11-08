#ifndef WORKERTHREADS_H_INCLUDED
#define WORKERTHREADS_H_INCLUDED

#include <exception>
#include <new>
#include <system_error>
#include <thread>
#include <vector>

#include "TGLException.h"

// Runs work() on the calling thread and on up to nthreads - 1 more threads, and returns once
// every one has finished. work() takes its own share of the job (an atomic counter of items,
// say), so it does not matter how many threads end up running it.
//
// work() must not call the R API nor read R's globals (NA_REAL, R_NaString, ...), which are
// not thread-safe: fetch whatever it needs from R on the calling thread first.
//
// A thread that cannot be started (std::system_error, for one under a process limit) leaves its
// share to the others; the calling thread always runs work(). An exception thrown by work() in
// any thread is caught there - one escaping a thread would end the process through
// std::terminate - and raised on the calling thread after all of them have finished: a
// TGLException or std::bad_alloc as it is, anything else as a TGLError with its message.
template <typename Work>
void run_in_threads(unsigned nthreads, const Work &work)
{
	if (!nthreads)
		nthreads = 1;

	std::vector<std::exception_ptr> errors(nthreads);
	auto guarded = [&](unsigned i) {
		try {
			work();
		} catch (...) {
			errors[i] = std::current_exception();
		}
	};

	std::vector<std::thread> threads;
	try {
		threads.reserve(nthreads - 1);
		for (unsigned i = 1; i < nthreads; ++i)
			threads.emplace_back(guarded, i);
	} catch (const std::system_error &) {
	} catch (const std::bad_alloc &) {
	}

	guarded(0);
	for (auto &thread : threads)
		thread.join();

	for (const std::exception_ptr &error : errors) {
		if (!error)
			continue;
		try {
			std::rethrow_exception(error);
		} catch (const TGLException &) {
			throw;
		} catch (const std::bad_alloc &) {
			throw;
		} catch (const std::exception &e) {
			TGLError("%s", e.what());
		} catch (...) {
			TGLError("A worker thread failed");
		}
	}
}

#endif
