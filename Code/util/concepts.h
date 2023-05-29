// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_CONCEPTS_H
#define HEMELB_UTIL_CONCEPTS_H

#include <concepts>
#include "util/traits.h"

namespace hemelb::util {

    template <typename T>
    concept optional = is_optional<T>::value;

    template<typename Base, typename Derived>
    concept base_of = std::derived_from<Derived, Base>;

    template<typename Base, typename Derived>
    concept pointer_base_of = std::derived_from<
            std::decay_t<decltype(*std::declval<Derived>())>,
            std::decay_t<decltype(*std::declval<Base>())>
    >;

}

#endif
