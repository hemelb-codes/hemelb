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
  void zeroOutFOld(geometry::LatticeData * const _latDat);
  // Initializes a lattice to empty distributions and forces, execpt at one
  // place
  template<class LATTICE>
    void allZeroButOne(
        geometry::LatticeData *_latDat, LatticeVector const &_pos);

  // FNew at given site
  template<class Lattice>
    distribn_t const * GetFNew(
        geometry::LatticeData *_latDat, LatticeVector const &_pos);

  // Linear distribution. Lambda pre-c++11.
  struct Linear {
    Linear(Dimensionless const _c[]) {
      for(size_t i(0); i < 4; ++i)
        c[i] = _c[i];
    }
    Linear(Dimensionless _o, Dimensionless _ax, Dimensionless _ay,
        Dimensionless _az) {
      c[0] = _o; c[1] = _ax; c[2] = _ay; c[3] = _az;
    }
    Linear(Dimensionless _o, LatticePosition const &_pos) {
      c[0] = _o; c[1] = _pos[0]; c[2] = _pos[1]; c[3] =_pos[2];
    }
    Dimensionless c[4];
    Dimensionless operator()(PhysicalVelocity const &_v) {
      return c[0] + _v[0] * c[1] + _v[1] * c[2] + _v[2] * c[3];
    }
  };

  // Population i set to a linear distribution
  // f_i(r) = _offset + (_ax, _ay, _az) * r
  template<class LATTICE>
    void setUpDistribution(
        geometry::LatticeData *_latDat, size_t _i, Linear const &_linear);
  template<class LATTICE>
    void setUpDistribution(
        geometry::LatticeData *_latDat, size_t _i, Dimensionless _offset,
        Dimensionless _ax, Dimensionless _ay, Dimensionless _az) {
      setUpDistribution<LATTICE>(_latDat, _i, Linear(_offset, _ax, _ay, _az));
    }


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

      // Sets up linear distribution f_i(coeffs) = (1, x, y, z) * coeffs
      // Where * is the dot product
      template<class LATTICE, class FUNCTION>
        void SetUpDistribution(size_t _i, FUNCTION function);

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
  inline void ZeroOutForces(geometry::LatticeData * const _latDat) {
    LatticeDataAccess(_latDat).ZeroOutForces();
  }

  template<class LATTICE, class FUNCTION>
    void LatticeDataAccess :: SetUpDistribution(size_t _i, FUNCTION function) {
      for(size_t i(0); i < latDat->GetLocalFluidSiteCount(); ++i) {
        geometry::Site<geometry::LatticeData> site = latDat->GetSite(i);
        LatticeVector const pos = site.GetGlobalSiteCoords();
        LatticePosition const pos_real(pos[0], pos[1], pos[2]);
        distribn_t const * const siteFOld(site.GetFOld<LATTICE>());
        distribn_t const * const firstFOld = &latDat->oldDistributions[0];
        size_t const indexFOld(siteFOld - firstFOld);
        latDat->oldDistributions[indexFOld + _i] = function(pos_real);
      }
    }


  template<class LATTICE>
    void allZeroButOne(
        geometry::LatticeData *_latDat, LatticeVector const &_pos) {
      LatticeDataAccess manip(_latDat);
      manip.ZeroOutFOld();
      manip.ZeroOutForces();
      for(size_t i(0); i < LATTICE::NUMVECTORS; ++i)
        manip.SetFOld<LATTICE>(_pos, i, 0.5 + i);
      _latDat->GetSite(_pos).SetForce(LatticeForceVector(1, 2, 3));
    }

  template<class LATTICE>
    distribn_t const * GetFNew(
        geometry::LatticeData *_latDat, LatticeVector const &_pos) {
      return LatticeDataAccess(_latDat).GetFNew<LATTICE>(_pos);
    }

  // Population i set to a linear distribution
  // f_i(r) = _offset + (_ax, _ay, _az) * r
  template<class LATTICE>
    void setUpDistribution(
        geometry::LatticeData *_latDat, size_t _i, Linear const &_linear) {
      LatticeDataAccess(_latDat)
        .SetUpDistribution<LATTICE, Linear>(_i, _linear);
    }


}}}
#endif
