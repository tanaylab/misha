#ifndef RUNNING_MAX_DEQUE_H
#define RUNNING_MAX_DEQUE_H

#include <deque>
#include <utility>
#include <limits>
#include <vector>

// Sliding-window maximum over floats with indices
struct RunningMaxDeque {
    // store pairs (index, value) with decreasing value
    std::deque<std::pair<size_t, float>> dq;
    size_t base_index = 0; // index of the first element in the window

    void clear() {
        dq.clear();
        base_index = 0;
    }

    void init(const std::vector<float>& vals) {
        clear();
        for (size_t k = 0; k < vals.size(); ++k) push(vals[k]);
    }

    // push value at next index (base_index + window_size_so_far)
    void push(float v) {
        size_t idx = base_index + size();
        while (!dq.empty() && dq.back().second <= v) dq.pop_back();
        dq.emplace_back(idx, v);
    }

    // pop one value from the front (advances base_index by 1)
    void pop_front() {
        if (!dq.empty() && dq.front().first == base_index) dq.pop_front();
        ++base_index;
    }

    size_t size() const {
        // window size = (last_idx - base_index + 1). We only need dq to eject front at base_index
        // but callers track size externally when needed.
        if (dq.empty()) return 0;
        return dq.back().first - base_index + 1; // not exact window size; callers shouldn't rely on this.
    }

    float value() const {
        if (dq.empty()) return -std::numeric_limits<float>::infinity();
        return dq.front().second;
    }

    // Argmax index within the logical stream (for MAX_POS)
    size_t argmax_index() const {
        return dq.empty() ? size_t(-1) : dq.front().first;
    }
};

#endif // RUNNING_MAX_DEQUE_H
