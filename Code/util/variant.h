// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_VARIANT_H
#define HEMELB_UTIL_VARIANT_H

#include <variant>

namespace hemelb
{
  namespace detail {
    // Helpers for allowing one to pass a number of lambdas to visit a
    // std:variant (see visit in namespaces below).
    template<typename... Ts>
    struct make_overload: Ts... {
      using Ts::operator()...;
    };

    template<typename... Ts>
    make_overload(Ts...) -> make_overload<Ts...>;
  }

  // Helper to visit a variant with a bunch of lamdbas. E.g.:
  //
  // overload_visit(src_var,
  //   [](first_type val) {/* whatever */},
  //   [](second_type val) { /* different thing for this type */},
  //   ... callables for other possible types ...
  // );
  template <typename Variant, typename... Callables>
  auto overload_visit(Variant&& v, Callables&&... calls) {
    return std::visit(detail::make_overload{std::forward<Callables>(calls)...}, v);
  }

}

#endif
