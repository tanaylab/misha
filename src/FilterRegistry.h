#ifndef FILTERREGISTRY_H_
#define FILTERREGISTRY_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <mutex>
#include "Filter.h"

/**
 * FilterRegistry: Singleton registry for compiled genomic filters.
 * Filters are cached by a unique key to avoid recompilation of identical masks.
 * Thread-safe for concurrent access.
 */
class FilterRegistry {
public:
    /**
     * Get the singleton instance.
     */
    static FilterRegistry& instance();

    /**
     * Get a filter by key. Returns nullptr if not found.
     */
    std::shared_ptr<Genome1DFilter> get(const std::string& key);

    /**
     * Store a filter under the given key.
     */
    void put(const std::string& key, std::shared_ptr<Genome1DFilter> filter);

    /**
     * Remove a filter by key. Returns true if found and removed.
     */
    bool remove(const std::string& key);

    /**
     * Clear all filters from the registry.
     */
    void clear();

    /**
     * Get number of registered filters.
     */
    size_t size() const;

private:
    FilterRegistry() = default;
    FilterRegistry(const FilterRegistry&) = delete;
    FilterRegistry& operator=(const FilterRegistry&) = delete;

    std::unordered_map<std::string, std::shared_ptr<Genome1DFilter>> filters_;
    mutable std::mutex mutex_;
};

#endif // FILTERREGISTRY_H_
