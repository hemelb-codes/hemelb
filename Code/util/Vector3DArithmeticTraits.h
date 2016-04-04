
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_VECTOR3DARITHMETICTRAITS_H
#define HEMELB_UTIL_VECTOR3DARITHMETICTRAITS_H

namespace hemelb
{
  namespace util
  {
    /**
     * This traits metastructure allows operator* and operator/ to choose the right
     * type for the result of a product or a division between arguments of potentially
     * different types.
     */
    template<class OperatorArgument1Type, class OperatorArgument2Type>
    struct Vector3DArithmeticTraits
    {
        /*
         * If this assertion trips, it means that Vector3D::operator* or Vector3D::operator/
         * have not been tested with this combination of types (OperatorArgument1Type, OperatorArgument2Type)
         *
         * One would like to write HEMELB_STATIC_ASSERT(false), however if the
         * static assertion is not dependent upon one or more template parameters, then
         * the compiler is permitted to evaluate the static assertion at the point it is
         * first seen, irrespective of whether the template is ever instantiated.
         *
         * See http://www.boost.org/doc/libs/1_49_0/doc/html/boost_staticassert.html
         */
        HEMELB_STATIC_ASSERT(sizeof(OperatorArgument1Type) == 0);

        /*
         * Boost alternative including an error message. If the C++0x static_assert feature is
         * not available, BOOST_STATIC_ASSERT_MSG(x, msg) will be treated as BOOST_STATIC_ASSERT(x)
         * and unfortunately you won't see the message.
         */
        //BOOST_STATIC_ASSERT_MSG(sizeof(T1) == 0, "Vector3D has not been tested with this combination of types");
    };

    // Trivial case: both arguments share type.
    template<class OperatorArgumentsType>
    struct Vector3DArithmeticTraits<OperatorArgumentsType, OperatorArgumentsType>
    {
        typedef OperatorArgumentsType operatorReturnType;
    };

    // The result must be stored as a floating point number.
    template<>
    struct Vector3DArithmeticTraits<int, float>
    {
        typedef float operatorReturnType;
    };

    // The result must be stored as a double-precision floating-point number.
    template<>
    struct Vector3DArithmeticTraits<int, double>
    {
        typedef double operatorReturnType;
    };

    // Keep the most precise type.
    template<>
    struct Vector3DArithmeticTraits<float, double>
    {
        typedef double operatorReturnType;
    };

    // This specialisation is needed by the setup tool. Keep the most precise type.
    template<>
    struct Vector3DArithmeticTraits<int, unsigned>
    {
        typedef unsigned operatorReturnType;
    };
    
    // This specialisation is needed by the setup tool. Keep the most precise type.
    template<>
    struct Vector3DArithmeticTraits<long int, int>
    {
        typedef long int operatorReturnType;
    };
    
    // This specialisation is needed by the setup tool. Keep the most precise type.
    template<>
    struct Vector3DArithmeticTraits<long long int, int>
    {
        typedef long long int operatorReturnType;
    };

  }
}
#endif
