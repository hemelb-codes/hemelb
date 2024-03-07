// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MPIDATATYPE_H
#define HEMELB_NET_MPIDATATYPE_H

#include <cstdint>
#include <cstddef>
#include <array>

#include <mpi.h>
#include <boost/hana/basic_tuple.hpp>
#include <boost/hana/map.hpp>
#include <boost/hana/pair.hpp>
#include <boost/hana/set.hpp>
#include <boost/hana/type.hpp>
#include <boost/hana/fwd/fold_left.hpp>
#include <boost/hana/at_key.hpp>

#include "net/MpiError.h"

// Implement C++ type to MPI_Datatype mapping
//
// User code (mainly in net::Mpi*) should just call MpiDataType<T>()
// to get the corresponding type.
//
// If the type requires creating and committing (with
// MPI_Type_commit), you should do this in a specialisation of
// MpiDataTypeRegistrationTraits. See below for an example which does
// this for std::array instances.
//
// If your type (T) can be treated as another (e.g. an enum), you can
// provide a function MpiDataType(T const&) which returns
// MpiDataType<Equivalent>()

namespace hemelb::net
{

    namespace detail {
        namespace hana = boost::hana;

        // The predefined types from the MPI standard, without
        // complex. We will use this twice below.
#define HEMELB_MPI_PREDEFINED_TYPELIST \
        HEMELB_MPI_BUILTIN_TYPE(char, MPI_CHAR), \
        HEMELB_MPI_BUILTIN_TYPE(short, MPI_SHORT), \
        HEMELB_MPI_BUILTIN_TYPE(int, MPI_INT), \
        HEMELB_MPI_BUILTIN_TYPE(long, MPI_LONG), \
        HEMELB_MPI_BUILTIN_TYPE(long long, MPI_LONG_LONG), \
        HEMELB_MPI_BUILTIN_TYPE(signed char, MPI_SIGNED_CHAR), \
        HEMELB_MPI_BUILTIN_TYPE(unsigned char, MPI_UNSIGNED_CHAR), \
        HEMELB_MPI_BUILTIN_TYPE(unsigned short, MPI_UNSIGNED_SHORT), \
        HEMELB_MPI_BUILTIN_TYPE(unsigned, MPI_UNSIGNED), \
        HEMELB_MPI_BUILTIN_TYPE(unsigned long, MPI_UNSIGNED_LONG), \
        HEMELB_MPI_BUILTIN_TYPE(unsigned long long, MPI_UNSIGNED_LONG_LONG), \
        HEMELB_MPI_BUILTIN_TYPE(float, MPI_FLOAT), \
        HEMELB_MPI_BUILTIN_TYPE(double, MPI_DOUBLE), \
        HEMELB_MPI_BUILTIN_TYPE(std::int8_t, MPI_INT8_T), \
        HEMELB_MPI_BUILTIN_TYPE(std::int16_t, MPI_INT16_T), \
        HEMELB_MPI_BUILTIN_TYPE(std::int32_t, MPI_INT32_T), \
        HEMELB_MPI_BUILTIN_TYPE(std::int64_t, MPI_INT64_T), \
        HEMELB_MPI_BUILTIN_TYPE(std::uint8_t, MPI_UINT8_T), \
        HEMELB_MPI_BUILTIN_TYPE(std::uint16_t, MPI_UINT16_T), \
        HEMELB_MPI_BUILTIN_TYPE(std::uint32_t, MPI_UINT32_T), \
        HEMELB_MPI_BUILTIN_TYPE(std::uint64_t, MPI_UINT64_T), \
        HEMELB_MPI_BUILTIN_TYPE(std::byte, MPI_BYTE)


        // We need a constexpr way to find out if a type is
        // predefined, whether or not we can get the MPI_Datatype that
        // it corresponds to
#define HEMELB_MPI_BUILTIN_TYPE(cpp, mpi) hana::type_c<cpp>

        constexpr auto MpiPredefinedTypes = hana::to_set(
            hana::make_basic_tuple(
                HEMELB_MPI_PREDEFINED_TYPELIST
            )
        );

#undef HEMELB_MPI_BUILTIN_TYPE


        // Accessing the datatype may or may not be constexpr, create
        // a map for this.

#define HEMELB_MPI_BUILTIN_TYPE(cpp, mpi) hana::make_pair(hana::type_c<cpp>, mpi)

        inline const auto MpiPredefinedTypeMap = hana::fold_left(
            hana::make_basic_tuple(
                HEMELB_MPI_PREDEFINED_TYPELIST
            ),
            hana::make_map(), hana::insert
        );

#undef HEMELB_MPI_BUILTIN_TYPE

        // Helper for below concept
        template <typename T>
        constexpr bool is_predefined() {
          constexpr auto ans = hana::contains(MpiPredefinedTypes, hana::type_c<T>);
          return decltype(ans)::value;
        }

    } // detail

    // Is the type one of those enumerated in the MPI standard as
    // being predefined?
    template <typename T>
    concept Predefined = detail::is_predefined<T>();

    // Primary template undefined - specialise if you need to commit
    template <typename T>
    struct MpiDataTypeRegistrationTraits;

    // Here we build the overload set to actually return the MPI type.
    // We provide 2 main overloads:
    // - for predefined types, use the machinery in detail
    // - for other types, use the RegistrationTraits above
    // You can provide overloads in other namespaces to implement
    // aliases (eg for enums)
    template<typename T>
    MPI_Datatype MpiDataType(T const&) {
        static MPI_Datatype DT = MPI_DATATYPE_NULL;
        if (DT == MPI_DATATYPE_NULL) {
            DT = MpiDataTypeRegistrationTraits<T>::Register();
        }
        return DT;
    }

    template<Predefined P>
    MPI_Datatype MpiDataType(P const&) {
        return detail::MpiPredefinedTypeMap[detail::hana::type_c<P>];
    }


    // User interface for getting the MPI_Datatype
    template<typename T>
    MPI_Datatype MpiDataType() {
        static_assert(
            std::is_trivial_v<T>,
            "MPI only works with trivially copyable types and we further require default construction"
        );
        return MpiDataType(T{});
    }

    // Template to register arrays
    template<typename T, std::size_t N>
    struct MpiDataTypeRegistrationTraits<std::array<T, N>> {
        static MPI_Datatype Register() {
            int blocklengths[1] = { N };
            MPI_Datatype types[1] = { MpiDataType<T>() };
            MPI_Aint displacements[1] = { 0 };
            MPI_Datatype ret;
            MpiCall{MPI_Type_create_struct}(1, blocklengths, displacements, types, &ret);
            MpiCall{MPI_Type_commit}(&ret);
            return ret;
        }
    };

}
#endif // HEMELB_NET_MPIDATATYPE_H
