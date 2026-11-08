#ifndef _SYSTEM_H_INCLUDED_
#define _SYSTEM_H_INCLUDED_
#include <cstdint>

// This function returns the actual number of unique bytes (cannot shared with other processes) used by pid.
// This is done by counting Private Dirty memory mappings.
// 
// Note: the result might significantly differ from what is being reported by "top" or "ps" utilities since they count
//       also the memory that can be shared between different processes.
// 
// Note: this function has high run-time cost since it needs to open and to parse system files.
// 
// !!!! COMPATIBILITY: only Linux, on other platforms returns 0. !!!!
uint64_t get_unique_mem_usage(pid_t pid = 0);

// Resident size of pid in bytes, from /proc/<pid>/statm: kernel counters, read without walking
// the page tables. It counts every resident page, shared ones included, so it is an upper
// bound of get_unique_mem_usage() (up to the few pages the kernel's counters may lag by).
//
// !!!! COMPATIBILITY: only Linux, on other platforms returns 0. !!!!
uint64_t get_resident_mem_usage(pid_t pid);

#endif
