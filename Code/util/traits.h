// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_TRAITS_H
#define HEMELB_UTIL_TRAITS_H

#include <type_traits>

namespace hemelb::util {

  // Type trait for a class being an instantiation of std::optional.
  template <typename T>
  struct is_optional : std::false_type {};
  template <typename T>
  struct is_optional<std::optional<T>> : std::true_type {};

  template <typename T>
  inline constexpr bool is_optional_v = is_optional<T>::value;

}

#endif
