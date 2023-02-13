// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_CT_STRING_H
#define HEMELB_UTIL_CT_STRING_H

#include <cstddef>
#include <string>
#include <type_traits>

namespace hemelb {

    // A compile-time string suitable for use as a non-type template parameter.
    // Easily comverts to C string or std::string.
    template<std::size_t N>
    struct ct_string {
        char str_[N + 1];

        constexpr ct_string(const char (&chars)[N + 1]) noexcept {
            for (std::size_t i = 0; i < N; ++i)
                str_[i] = chars[i];
            str_[N] = 0;
        }

        constexpr std::string str() const {
            return {str_, N};
        }

        constexpr char const *c_str() const {
            return str_;
        }

        template<std::size_t M>
        friend constexpr auto operator<=>(const ct_string &left,
                                          const ct_string<M> &right) {
            return left.str().compare(right.str()) <=> 0;
        }
    };

    template<std::size_t N>
    ct_string(const char (&str)[N]) -> ct_string<N - 1>;

    template<typename T>
    struct is_ct_string : std::false_type {
    };
    template<std::size_t N>
    struct is_ct_string<ct_string<N>> : std::true_type {
    };
    template<typename T>
    concept is_ct_string_v = is_ct_string<std::decay_t<T>>::value;

}
#endif