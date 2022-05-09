// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELBSETUPTOOL_RANGE_HPP
#define HEMELBSETUPTOOL_RANGE_HPP

#include <concepts>
#include <ranges>

// Simple use like python's `range` generator
//
// for (auto i: range(5)) std::cout << i << ",";
//
// Will print 0,1,2,3,4,
// 
template <std::integral T>
constexpr auto range(T lo, T hi) {
  return std::ranges::iota_view<T,T>{lo, hi};
}
template <std::integral T>
constexpr auto range(T hi) {
  return std::ranges::iota_view<T,T>{0, hi};
}

#endif
