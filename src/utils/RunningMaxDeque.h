#ifndef RUNNING_MAX_DEQUE_H
#define RUNNING_MAX_DEQUE_H

#include <deque>
#include <utility>
#include <limits>
#include <vector>

// Sliding-window maximum over floats with genomic positions and strand direction
struct RunningMaxDeque {
    // store triples (genomic_position, value, direction) with decreasing value
    struct Entry {
        int64_t genomic_pos;  // Absolute genomic coordinate (not target-relative)
        float val;
        int dir;  // 1 for forward strand, -1 for reverse strand

        Entry(int64_t pos, float v, int d) : genomic_pos(pos), val(v), dir(d) {}
    };

    std::deque<Entry> dq;
    int64_t base_genomic_pos = -1; // genomic position of the first element in the window

    void clear() {
        dq.clear();
        base_genomic_pos = -1;
    }

    void init(const std::vector<float>& vals, const std::vector<int>& dirs, const std::vector<int64_t>& genomic_positions) {
        clear();
        for (size_t k = 0; k < vals.size(); ++k) {
            int dir = (k < dirs.size()) ? dirs[k] : 1;
            int64_t pos = (k < genomic_positions.size()) ? genomic_positions[k] : int64_t(k);
            push(vals[k], dir, pos);
        }
    }

    void init(const std::vector<float>& vals, const std::vector<int>& dirs) {
        // For backward compatibility - use sequential indices as positions
        std::vector<int64_t> positions(vals.size());
        for (size_t i = 0; i < vals.size(); ++i) positions[i] = int64_t(i);
        init(vals, dirs, positions);
    }

    void init(const std::vector<float>& vals) {
        // For backward compatibility - assume forward strand, sequential positions
        std::vector<int> dirs(vals.size(), 1);
        std::vector<int64_t> positions(vals.size());
        for (size_t i = 0; i < vals.size(); ++i) positions[i] = int64_t(i);
        init(vals, dirs, positions);
    }

    // push value with genomic position and direction
    void push(float v, int dir, int64_t genomic_pos) {
        if (dq.empty()) {
            base_genomic_pos = genomic_pos;
        }
        while (!dq.empty() && dq.back().val < v) dq.pop_back();
        dq.emplace_back(genomic_pos, v, dir);
    }

    // Backward compatibility: push without explicit genomic position
    void push(float v, int dir = 1) {
        int64_t pos = dq.empty() ? 0 : (dq.back().genomic_pos + 1);
        push(v, dir, pos);
    }

    // pop all entries before new_base_genomic_pos (for sliding window)
    void pop_front(int64_t new_base_genomic_pos) {
        while (!dq.empty() && dq.front().genomic_pos < new_base_genomic_pos) {
            dq.pop_front();
        }
        base_genomic_pos = new_base_genomic_pos;
    }

    // Backward compatibility: pop one position
    void pop_front() {
        pop_front(base_genomic_pos + 1);
    }

    size_t size() const {
        // Return number of elements in deque (not necessarily window size)
        return dq.size();
    }

    float value() const {
        if (dq.empty()) return -std::numeric_limits<float>::infinity();
        return dq.front().val;
    }

    // Argmax genomic position (for MAX_POS)
    int64_t argmax_genomic_position() const {
        return dq.empty() ? -1 : dq.front().genomic_pos;
    }

    // Backward compatibility: return as size_t (deprecated)
    size_t argmax_index() const {
        int64_t pos = argmax_genomic_position();
        return (pos < 0) ? size_t(-1) : size_t(pos);
    }

    // Get the direction of the argmax (for MAX_POS with bidirectional PSSMs)
    int argmax_direction() const {
        return dq.empty() ? 1 : dq.front().dir;
    }
};

#endif // RUNNING_MAX_DEQUE_H
