// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_LATTICES_H
#define HEMELB_LB_LATTICES_H

#include "lb/lattices/D3Q15.h"
#include "lb/lattices/D3Q19.h"
#include "lb/lattices/D3Q27.h"
#include "lb/lattices/D3Q15i.h"
#include "build_info.h"

namespace hemelb::lb {

    namespace detail {
        constexpr auto get_default_lattice() {
            constexpr auto LAT = build_info::LATTICE;
            if constexpr (LAT == "D3Q15") {
                return D3Q15{};
            } else if constexpr (LAT == "D3Q19") {
                return D3Q19{};
            } else if constexpr (LAT == "D3Q27") {
                return D3Q27{};
            } else if constexpr (LAT == "D3Q15i") {
                return D3Q15i{};
            } else {
                throw (Exception() << "Configured with invalid LATTICE");
            }
        }
    }
    using DefaultLattice = decltype(detail::get_default_lattice());
}

#endif /* HEMELB_LB_LATTICES_H */
