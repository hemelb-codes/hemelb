// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_MOMENTBASIS_BASIS_HELPERS_H
#define HEMELB_LB_KERNELS_MOMENTBASIS_BASIS_HELPERS_H

#include <array>
#include <numeric>

namespace hemelb::lb::kernels::momentBasis {

    template<typename T, std::size_t NUM_MOMS, std::size_t NUM_VELS>
    struct moment_traits {
        using MatrixType = std::array<std::array<T, NUM_VELS>, NUM_MOMS>;

        static constexpr auto DiagSelfProduct(MatrixType const &mat) {
            std::array<T, NUM_MOMS> ans;
            for (unsigned i = 0; i < NUM_MOMS; ++i) {
                auto &row = mat[i];
                ans[i] = std::inner_product(row.begin(), row.end(), row.begin(), T{0});
            }
            return ans;
        }
    };

}
#endif