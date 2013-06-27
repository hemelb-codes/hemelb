//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UTIL_BESSEL_H
#define HEMELB_UTIL_BESSEL_H

#include <complex>

namespace hemelb
{
  namespace util
  {

    /**
     * Evaluates the Bessel function of the first kind order 0 at z. The solution is approximated
     * with numTerms terms of the series
     *
     *    J0(z) ~=  1 - 0.25*z*z/1!^2 + (0.25*z*z)^2/2!^2 - (0.25*z*z)^3/3!^2 + ...
     *
     * See formula 9.1.12 of Abramowitz and Stegun, Handbook of Mathematical Functions.
     *
     * At the time of writing the Boost implementation of the Bessel functions did not support
     * complex arguments.
     *
     * @param z point in the complex plane where we want to evaluate J0
     * @param numTerms number of terms in the series expansion used to approximate the function
     * @return J0(z)
     */
    std::complex<double>
    BesselJ0ComplexArgument(const std::complex<double>& z, double tolSq = 1e-12);
  }
}

#endif // HEMELB_UTIL_BESSEL_H
