// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_SPAN_H
#define HEMELB_UTIL_SPAN_H

#include <span>
#include <vector>

namespace hemelb {
    // Helpers to turn vectors into spans
    template <typename T>
    std::span<T> to_span(std::vector<T>& v) {
        return {v.data(), v.size()};
    }

    template <typename T>
    std::span<T const> to_const_span(std::vector<T> const& v) {
        return {v.data(), v.size()};
    }

    template <typename T>
    std::span<T const> to_span(std::vector<T> const& v) {
        return {v.data(), v.size()};
    }
}

#endif
