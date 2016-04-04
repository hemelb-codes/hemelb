
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_BESSEL_H
#define HEMELB_UTIL_BESSEL_H

#include <complex>

namespace hemelb
{
  namespace util
  {

    /**
     * Evaluates the Bessel function of the first kind order 0 at z. The solution is approximated
     * by the series
     *
     *    J0(z) ~=  1 - 0.25*z*z/1!^2 + (0.25*z*z)^2/2!^2 - (0.25*z*z)^3/3!^2 + ...
     *
     * with the termination condition that the norm of the last term is less than
     * the specified tolerance.
     *
     * See formula 9.1.12 of Abramowitz and Stegun, Handbook of Mathematical Functions.
     *
     * At the time of writing the Boost implementation of the Bessel functions did not support
     * complex arguments.
     *
     * @param z point in the complex plane where we want to evaluate J0
     * @param tolSq the square of the tolerance used to terminate the expansion.
     * @return J0(z)
     */
    std::complex<double>
    BesselJ0ComplexArgument(const std::complex<double>& z, double tolSq = 1e-12);
  }
}

#endif // HEMELB_UTIL_BESSEL_H
