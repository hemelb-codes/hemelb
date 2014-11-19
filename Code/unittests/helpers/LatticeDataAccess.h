//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_HELPERS_LATTICEDATAACCESS_H
#define HEMELB_UNITTESTS_HELPERS_LATTICEDATAACCESS_H

#include <algorithm>
#include "geometry/LatticeData.h"

namespace hemelb { namespace unittests { namespace helpers {

  // Empty FOld distribution
  void ZeroOutFOld(geometry::LatticeData * const _latDat);

  // Class to setup/manipulate lattice data
  class LatticeDataAccess {
    public:
      LatticeDataAccess(geometry::LatticeData * const _latDat)
        : latDat(_latDat) {};

      // Empty distribution
      void ZeroOutFOld() const;
      // Empty Forces
      void ZeroOutForces() const;
      // Set FOld for a given site and direction
      template<class LATTICE>
        void SetFOld(
            site_t _x, site_t _y, site_t _z, site_t _dir, distribn_t _value
        ) const {
          SetFOld<LATTICE>(LatticeVector(_x, _y, _z), _dir, _value);
        }
      template<class LATTICE>
        void SetFOld(
            LatticeVector const &_pos, site_t _dir, distribn_t _value) const;

      // Get FNew for a given site and direction
      template<class LATTICE>
        distribn_t const * GetFNew(site_t _x, site_t _y, site_t _z) const {
          return GetFNew<LATTICE>(LatticeVector(_x, _y, _z));
        }
      template<class LATTICE>
        distribn_t const * GetFNew(LatticeVector const &_pos) const;

    protected:
      // Reference to the wrapped lattice data.
      geometry::LatticeData * const latDat;
  };

  inline void LatticeDataAccess :: ZeroOutFOld() const {
    std::fill(
        latDat->oldDistributions.begin(), latDat->oldDistributions.end(), 0);
  }
  inline void LatticeDataAccess :: ZeroOutForces() const {
    std::fill(latDat->forceAtSite.begin(), latDat->forceAtSite.end(), 0);
  }

  template<class LATTICE>
    void LatticeDataAccess :: SetFOld(
        LatticeVector const &_pos, site_t _dir, distribn_t _value) const {
      // Figure location of distribution in memory using const access to FOld.
      // This way, access is resilient versus (some) changes in memory layout.
      geometry::Site<geometry::LatticeData> const site(latDat->GetSite(_pos));
      distribn_t const * const siteFOld(site.GetFOld<LATTICE>());
      distribn_t const * const firstFOld = &latDat->oldDistributions[0];
      size_t const indexFOld(siteFOld - firstFOld);
      latDat->oldDistributions[indexFOld + _dir] = _value;
    }

  template<class LATTICE>
    distribn_t const *
    LatticeDataAccess :: GetFNew(LatticeVector const &_pos) const {
      // Figure location of distribution in memory using const access to FOld.
      // This way, access is resilient versus (some) changes in memory layout.
      geometry::Site<geometry::LatticeData> const site(latDat->GetSite(_pos));
      distribn_t const * const siteFOld(site.GetFOld<LATTICE>());
      distribn_t const * const firstFOld = &latDat->oldDistributions[0];
      size_t const indexFOld(siteFOld - firstFOld);
      return latDat->GetFNew(indexFOld);
    }
  inline void ZeroOutFOld(geometry::LatticeData * const _latDat) {
    LatticeDataAccess(_latDat).ZeroOutFOld();
  }

}}}
#endif
