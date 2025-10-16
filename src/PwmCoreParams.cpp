#include "PwmCoreParams.h"

#include <algorithm>
#include <cmath>

void PwmCoreParams::set_spatial_factors(const std::vector<float>& factors, int bin_size) {
    spat_factor = factors;
    spat_bin_size = std::max(1, bin_size);

    spat_log_factors.resize(spat_factor.size());
    for (size_t i = 0; i < spat_factor.size(); ++i) {
        float val = spat_factor[i];
        if (val <= 0.0f) {
            val = 1e-30f;
        }
        spat_log_factors[i] = std::log(val);
    }
}

void PwmCoreParams::apply_prior() {
    if (prior <= 0.0)
        return;
    pssm.add_dirichlet_prior(static_cast<float>(prior));
}
