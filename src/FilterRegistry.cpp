#include "FilterRegistry.h"
#include <algorithm>

FilterRegistry& FilterRegistry::instance() {
    static FilterRegistry instance;
    return instance;
}

std::shared_ptr<Genome1DFilter> FilterRegistry::get(const std::string& key) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = filters_.find(key);
    if (it != filters_.end()) {
        return it->second;
    }
    return nullptr;
}

void FilterRegistry::put(const std::string& key, std::shared_ptr<Genome1DFilter> filter) {
    std::lock_guard<std::mutex> lock(mutex_);
    filters_[key] = filter;
}

bool FilterRegistry::remove(const std::string& key) {
    std::lock_guard<std::mutex> lock(mutex_);
    return filters_.erase(key) > 0;
}

void FilterRegistry::clear() {
    std::lock_guard<std::mutex> lock(mutex_);
    filters_.clear();
}

size_t FilterRegistry::size() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return filters_.size();
}
