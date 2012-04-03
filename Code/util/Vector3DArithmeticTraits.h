#ifndef HEMELB_UTIL_VECTOR3DARITHMETICTRAITS_H
#define HEMELB_UTIL_VECTOR3DARITHMETICTRAITS_H

//#include <boost/static_assert.hpp>
#include "units.h"

namespace hemelb
{
  namespace util
  {
    /**
     * This traits metastructure allows operator* to choose the right type for the
     * result of a product
     */
    template<class T1, class T2>
    struct Vector3DArithmeticTraits
    {
        /*
         * One would like to write BOOST_STATIC_ASSERT_MSG(false,"..."), however if the
         * static assertion is not dependent upon one or more template parameters, then
         * the compiler is permitted to evaluate the static assertion at the point it is
         * first seen, irrespective of whether the template is ever instantiated.
         *
         * If the C++0x static_assert feature is not available, BOOST_STATIC_ASSERT_MSG(x, msg)
         * will be treated as BOOST_STATIC_ASSERT(x) and unfortunately you won't see the message
         *
         * See http://www.boost.org/doc/libs/1_49_0/doc/html/boost_staticassert.html
         */

        ///@todo 138, reenable once BOOST is usable in other machines
        //BOOST_STATIC_ASSERT_MSG(sizeof(T1) == 0, "Vector3D has not been tested with this combination of types");
        HEMELB_STATIC_ASSERT(sizeof(T1) == 0);
        ///@todo 138, remove once BOOST is usable in other machines
        //typedef T1 type;
    };

    // Trivial case: both arguments share type.
    template<class T>
    struct Vector3DArithmeticTraits<T, T>
    {
        typedef T type;
    };

    // The result must be stored as a floating point number.
    template<>
    struct Vector3DArithmeticTraits<int, float>
    {
        typedef float type;
    };

    // Keep the most precise type.
    template<>
    struct Vector3DArithmeticTraits<float, double>
    {
        typedef double type;
    };

    // This specialisation is needed by the setup tool. Keep the most precise type.
    template<>
    struct Vector3DArithmeticTraits<int, unsigned>
    {
        typedef unsigned type;
    };

  }
}
#endif
