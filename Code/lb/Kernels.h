// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_H
#define HEMELB_LB_KERNELS_H

#include "build_info.h"

#include "lb/kernels/EntropicAnsumali.h"
#include "lb/kernels/EntropicChik.h"
#include "lb/kernels/LBGK.h"
#include "lb/kernels/LBGKNN.h"
#include "lb/kernels/MRT.h"
#include "lb/kernels/GuoForcingLBGK.h"
#include "lb/kernels/TRT.h"
#include "lb/kernels/MomentBases.h"
#include "lb/kernels/RheologyModels.h"

namespace hemelb::lb {
    namespace detail {
        template <lattice_type L>
        constexpr auto get_default_kernel(InitParams& i) {
            constexpr auto KERN = build_info::KERNEL;
            if constexpr (KERN == "LBGK") {
                return LBGK<L>{i};
            } else if constexpr (KERN == "EntropicAnsumali") {
                return EntropicAnsumali<L>{i};
            } else if constexpr (KERN == "EntropicChik") {
                return EntropicChik<L>{i};
            } else if constexpr (KERN == "MRT") {
                if constexpr (std::same_as<L, D3Q15>) {
                    return MRT<DHumieresD3Q15MRTBasis>{i};
                } else if constexpr(std::same_as<L, D3Q19>) {
                    return MRT<DHumieresD3Q19MRTBasis>{i};
                } else {
                    throw (Exception() << "No MRT basis for lattice");
                }
            } else if constexpr (KERN == "TRT") {
                return TRT<L>{i};
            } else if constexpr (KERN == "NNCY") {
                return LBGKNN<CarreauYasudaRheologyModelHumanFit, L>{i};
            } else if constexpr (KERN == "NNCYMOUSE") {
                return LBGKNN<CarreauYasudaRheologyModelMouseFit, L>{i};
            } else if constexpr (KERN == "NNC") {
                return LBGKNN<CassonRheologyModel, L>{i};
            } else if constexpr (KERN == "NNTPL") {
                return LBGKNN<TruncatedPowerLawRheologyModel, L>{i};
            } else if constexpr (KERN == "GuoForcingLBGK") {
                return GuoForcingLBGK<L>{i};
            } else {
                throw (Exception() << "Configured with invalid KERNEL");
            }
        }
    }

    template <lattice_type L>
    using DefaultKernel = decltype(detail::get_default_kernel<L>(std::declval<InitParams&>()));
}
#endif /* HEMELB_LB_KERNELS_H */
