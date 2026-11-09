#include <cstdint>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "BufferedFile.h"
#include "System.h"
#include "port.h"
#include "strutil.h"

uint64_t get_unique_mem_usage(pid_t pid)
{
	uint64_t mem_usage = 0;
	BufferedFile bf;
	char filename[100];
	vector<string> fields;

	snprintf(filename, sizeof(filename), "/proc/%ld/smaps", (long)pid);

	// count only private dirty bytes under heap section
	if (!bf.open(filename, "r")) {
		bool is_heap = false;

		while (1) {
			split_line_by_space_chars(bf, fields, 2);

			if (fields.empty()) 
				break;

			if (fields.size() == 6) {
				if (is_heap) {
					if (!fields[5].empty() && fields[5] != "[heap]")
						break;
				} else if (fields[5] == "[heap]") 
					is_heap = true;
			} else if (is_heap && fields.size() == 3 && fields[0] == "Private_Dirty:") {
				char *endptr;
				long num;

				num = strtol(fields[1].c_str(), &endptr, 10);
				if (!*endptr)
					mem_usage += num;
			}
		}
	}
	return mem_usage * 1024;
}

uint64_t get_resident_mem_usage(pid_t pid)
{
	char filename[100];
	snprintf(filename, sizeof(filename), "/proc/%ld/statm", (long)pid);

	FILE *fp = fopen(filename, "r");
	if (!fp)
		return 0;

	unsigned long size = 0, resident = 0;
	int n = fscanf(fp, "%lu %lu", &size, &resident);
	fclose(fp);

	long page_size = sysconf(_SC_PAGESIZE);
	if (n != 2 || page_size <= 0)
		return 0;
	return (uint64_t)resident * (uint64_t)page_size;
}

uint64_t get_anon_resident_mem_usage(pid_t pid)
{
	char filename[100];
	snprintf(filename, sizeof(filename), "/proc/%ld/statm", (long)pid);

	FILE *fp = fopen(filename, "r");
	if (!fp)
		return 0;

	unsigned long size = 0, resident = 0, shared = 0;
	int n = fscanf(fp, "%lu %lu %lu", &size, &resident, &shared);
	fclose(fp);

	long page_size = sysconf(_SC_PAGESIZE);
	if (n != 3 || page_size <= 0 || resident < shared)
		return 0;
	return (uint64_t)(resident - shared) * (uint64_t)page_size;
}

uint64_t get_anon_huge_pages(pid_t pid)
{
	char filename[100];
	char line[256];
	uint64_t kb = 0;

	snprintf(filename, sizeof(filename), "/proc/%ld/smaps_rollup", (long)pid);
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		snprintf(filename, sizeof(filename), "/proc/%ld/smaps", (long)pid);
		fp = fopen(filename, "r");
		if (!fp)
			return 0;
	}

	while (fgets(line, sizeof(line), fp)) {
		if (!strncmp(line, "AnonHugePages:", 14))
			kb += strtoull(line + 14, NULL, 10);
	}
	fclose(fp);
	return kb * 1024;
}
