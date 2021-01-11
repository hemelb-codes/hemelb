// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_COMPARISONS_H
#define HEMELB_UNITTESTS_HELPERS_COMPARISONS_H
namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {

      inline bool is_zero(util::Vector3D<double> const &_in, double _tol = 1e-8)
      {
        return std::sqrt(_in.GetMagnitudeSquared()) < _tol;
      }
      inline bool is_zero(double const _in, double _tol = 1e-8)
      {
        return std::abs(_in) < _tol;
      }

    }
  }
}
#endif
