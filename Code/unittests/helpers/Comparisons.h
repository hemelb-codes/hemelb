//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_COMPARISONS_H
#define HEMELB_UNITTESTS_COMPARISONS_H
namespace hemelb { namespace unittests { namespace helpers {

  inline bool is_zero(util::Vector3D<double> const &_in, double _tol = 1e-8) {
    return std::sqrt(_in.GetMagnitudeSquared()) < _tol;
  }
  inline  bool is_zero(double const _in, double _tol = 1e-8) {
    return std::abs(_in) < _tol;
  }

}}}
#endif
