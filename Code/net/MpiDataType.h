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
#include "net/MpiError.h"

namespace hemelb::net
{
    /*
     * There are templated functions for getting the MPI_Datatype
     * corresponding to a type or variable below. Use as:
     *     Foo bar;
     *     MPI_Datatype tp = MpiDataType(bar);
     *     // or alternatively
     *     tp = MpiDataType<Foo>();
     *
     * The generic version, for custom classes, is implemented using traits.
     * You must specialize the template, e.g.
     *
     *     template<>
     *     MPI_Datatype hemelb::net::MpiDataTypeTraits<Foo>::RegisterMpiDataType()
     *     {
     *       MPI_Datatype type;
     *       // Create the type
     *       MpiCall{MPI_Type_commit}(&type);
     *       return type;
     *     }
     *
     * Built-in MPI types have specializations defined which are also constexpr.
     * Below is a partial specialisation for std::array<T, N>.
     *
     * Primary template
     */
    template<typename T>
    class MpiDataTypeTraits
    {
      public:
        static const MPI_Datatype& GetMpiDataType()
        {
	  static_assert(!std::same_as<T, std::vector<unsigned>>);
          if (mpiType == MPI_DATATYPE_NULL)
          {
            mpiType = RegisterMpiDataType();
          }
          return mpiType;
        }
      private:
        static MPI_Datatype mpiType;
        static MPI_Datatype RegisterMpiDataType();
        template<class U> friend class MpiDataTypeTraits;
    };

    // Define the initial value
    template<typename T>
    MPI_Datatype MpiDataTypeTraits<T>::mpiType = MPI_DATATYPE_NULL;

    // Generic getters of data type info
    template<typename T>
    auto&& MpiDataType()
    {
      return MpiDataTypeTraits<T>::GetMpiDataType();
    }
    template<typename T>
    auto&& MpiDataType(const T&)
    {
      return MpiDataTypeTraits<T>::GetMpiDataType();
    }

    // Declare a constexpr specialisation for MPI built in types.
    #define HEMELB_MPI_BUILTIN_TYPE(T, MPI)                     \
    template<>                                                  \
    class MpiDataTypeTraits<T> {                                \
    public:                                                     \
        static constexpr MPI_Datatype const& GetMpiDataType() { \
            return mpiType;                                     \
        }                                                       \
    private:                                                    \
        static constexpr MPI_Datatype mpiType = MPI;            \
        template<class U> friend class MpiDataTypeTraits;       \
    }

    HEMELB_MPI_BUILTIN_TYPE(char, MPI_CHAR);
    HEMELB_MPI_BUILTIN_TYPE(signed char, MPI_SIGNED_CHAR);
    HEMELB_MPI_BUILTIN_TYPE(std::int16_t, MPI_SHORT);
    HEMELB_MPI_BUILTIN_TYPE(std::int32_t, MPI_INT);
    HEMELB_MPI_BUILTIN_TYPE(std::int64_t, MPI_LONG_LONG);
    // HEMELB_MPI_BUILTIN_TYPE(std::size_t, MPI_UNSIGNED_LONG);
    HEMELB_MPI_BUILTIN_TYPE(unsigned char, MPI_UNSIGNED_CHAR);
    HEMELB_MPI_BUILTIN_TYPE(std::uint16_t, MPI_UNSIGNED_SHORT);
    HEMELB_MPI_BUILTIN_TYPE(std::uint32_t, MPI_UNSIGNED);
    HEMELB_MPI_BUILTIN_TYPE(std::uint64_t, MPI_UNSIGNED_LONG_LONG);
    HEMELB_MPI_BUILTIN_TYPE(float, MPI_FLOAT);
    HEMELB_MPI_BUILTIN_TYPE(double, MPI_DOUBLE);

    // Specialisation for statically size arrays
    template <typename T, std::size_t N>
    class MpiDataTypeTraits<std::array<T, N>>
    {
    public:
        static const MPI_Datatype& GetMpiDataType() {
          if (mpiType == MPI_DATATYPE_NULL)
          {
            mpiType = RegisterMpiDataType();
          }
          return mpiType;
        }
    private:
        static MPI_Datatype mpiType;
        static MPI_Datatype RegisterMpiDataType() {
	  int blocklengths[1] = { N };
	  MPI_Datatype types[1] = { MpiDataType<T>() };
	  MPI_Aint displacements[1] = { 0 };
	  MPI_Datatype ret;
	  MpiCall{MPI_Type_create_struct}(1, blocklengths, displacements, types, &ret);
	  MpiCall{MPI_Type_commit}(&ret);
	  return ret;
	}

        template<class U> friend class MpiDataTypeTraits;
    };

    template <typename T, std::size_t N>
    MPI_Datatype MpiDataTypeTraits<std::array<T, N>>::mpiType = MPI_DATATYPE_NULL;

}
#endif // HEMELB_NET_MPIDATATYPE_H
