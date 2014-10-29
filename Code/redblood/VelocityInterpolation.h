//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_VELOCITY_INTERPOLATION_H
#define HEMELB_REDBLOOD_VELOCITY_INTERPOLATION_H

#include "units.h"
#include "Exception.h"
#include "redblood/Interpolation.h"
#include "redblood/stencil.h"

#include <vector>
#include <cassert>
#include <boost/shared_array.hpp>

namespace hemelb { namespace redblood {

  //! \brief Returns the interpolated velocity at given point
  //! \param[in] _center: off lattice position for which to interpolate.
  //! \param[in] _latdat: lattice data
  //! \param[in] _stencil: stencil to use for interpolation
  //! \tparam KERNEL: Needed to compute velocity
  template<class KERNEL>
    PhysicalVelocity interpolateVelocity(PhysicalPosition const &_center,
        geometry::LatticeData const &_latdat,
        stencil::types _stencil = stencil::FOUR_POINT);

  namespace details {
    // Computes velocity for a given index on the lattice
    template<class KERNEL> struct VelocityFromLatticeData {
      VelocityFromLatticeData(geometry::LatticeData const &_latdata)
            : latticeData(_latdata) {}
      // Computes velocity for given KERNEL and LatticeData
      LatticeVelocity operator()(LatticeVector const &_indices) const {
        typedef typename KERNEL::KHydroVars HydroVars;
        typedef typename KERNEL::LatticeType LatticeType;
        site_t const index(latticeData.GetContiguousSiteId(_indices));
        LatticeVelocity result;
        distribn_t density;
        LatticeType::CalculateDensityAndMomentum(
            latticeData.GetSite(index).template GetFOld<LatticeType>(),
            density, result.x, result.y, result.z
        );
        return result;
      }
      protected:
        //! Lattice data that holds the grid of population and forces
        geometry::LatticeData const & latticeData;
    };
  }

  template<class KERNEL>
    PhysicalVelocity interpolateVelocity(
        geometry::LatticeData const &_latdat, PhysicalPosition const &_center,
        stencil::types _stencil = stencil::FOUR_POINT) {
      return interpolate(
          details::VelocityFromLatticeData<KERNEL>(_latdat),
          interpolationIterator(_center, _stencil)
      );
    }
}} // hemelb::redblood
#endif
