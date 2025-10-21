#ifndef RUNNING_LOG_SUM_EXP_H
#define RUNNING_LOG_SUM_EXP_H

#include <deque>
#include <cmath>
#include <limits>
#include <vector>

// Maintains logsumexp over a fixed-size sliding window in O(1) amortized time.
// Numerically stable via max-trick and single rescaling when max changes.
struct RunningLogSumExp {
    double M = -std::numeric_limits<double>::infinity(); // current max
    double sum_scaled = 0.0;                              // sum exp(x - M)
    std::deque<float> window;                             // raw values (ℓ_i)
    std::deque<float> maxdq;                              // monotonic decreasing queue of candidates for M
    size_t W = 0;                                         // target window size (optional)

    void clear() {
        M = -std::numeric_limits<double>::infinity();
        sum_scaled = 0.0;
        window.clear();
        maxdq.clear();
        W = 0;
    }

    // Initialize from a full window of values.
    void init(const std::vector<float>& vals) {
        clear();
        if (vals.empty()) return;
        W = vals.size();
        // Find M
        for (float v : vals) if (v > M) M = v;
        // Build sum_scaled and maxdq
        for (float v : vals) {
            sum_scaled += std::exp(double(v) - M);
            while (!maxdq.empty() && maxdq.back() < v) maxdq.pop_back();
            maxdq.push_back(v);
            window.push_back(v);
        }
    }

    inline void push(float x) {
        // If we track W and caller wants fixed size behavior, they should pop beforehand.
        if (x > M) {
            // Rescale once
            if (std::isfinite(M)) sum_scaled *= std::exp(M - double(x));
            M = x;
        }
        sum_scaled += std::exp(double(x) - M);
        while (!maxdq.empty() && maxdq.back() < x) maxdq.pop_back();
        maxdq.push_back(x);
        window.push_back(x);
    }

    inline void pop_front() {
        float x = window.front();
        window.pop_front();
        // remove x from maxdq front if it matches
        if (!maxdq.empty() && maxdq.front() == x) {
            maxdq.pop_front();
            // If max changed, rescale to the new max once
            double M_new = maxdq.empty() ? -std::numeric_limits<double>::infinity()
                                         : double(maxdq.front());
            if (M_new != M) {
                if (std::isfinite(M_new)) {
                    sum_scaled *= std::exp(M - M_new);
                } else {
                    sum_scaled = 0.0;
                }
                M = M_new;
            }
        }
        // subtract x's scaled contribution (safe even if M changed due to rescale above)
        if (std::isfinite(M)) sum_scaled -= std::exp(double(x) - M);
        else sum_scaled = 0.0;
        if (sum_scaled < 0) sum_scaled = 0; // guard tiny negatives due to FP
    }

    inline double value() const {
        if (!std::isfinite(M)) return -std::numeric_limits<double>::infinity();
        return M + std::log(sum_scaled);
    }
};

#endif // RUNNING_LOG_SUM_EXP_H
