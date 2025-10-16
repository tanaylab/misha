#ifndef PWM_CORE_PARAMS_H_
#define PWM_CORE_PARAMS_H_

#include <vector>

#include "port.h"
#include "DnaPSSM.h"

struct PwmCoreParams {
    DnaPSSM pssm;
    bool bidirect = true;
    int strand_mode = 0;          // -1, 0, or 1
    double score_thresh = 0.0;
    int extend = 0;               // allowed extension of motif start positions (bases)
    std::vector<float> spat_factor;
    std::vector<float> spat_log_factors;
    int spat_bin_size = 1;
    double prior = 0.0;

    void set_spatial_factors(const std::vector<float>& factors, int bin_size);
    void apply_prior();
};

#endif // PWM_CORE_PARAMS_H_
