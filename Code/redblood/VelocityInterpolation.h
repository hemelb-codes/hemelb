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
    PhysicalVelocity interpolateVelocity(PhysicalPosition const &_center
        LatticeData const &_latdat,
        stencil::type _stencil = stencil::FOUR_POINT);

  namespace details {
    // Computes velocity for a given index on the lattice
    template<class KERNEL> struct VelocityFromLatticeData {
      Velocity(LatticeData const &_latdata) : latticeData(_latdat) {}
      // Computes velocity for given KERNEL and LatticeData
      LatticeVelocity operator(LatticeVector const &_indices) const {
        typedef kernels::HydroVars<KERNEL> HydroVars;
        site_t const index(latticeData.GetContiguousSiteId(_indices));
        HydroVars hydroVars(latticeData.GetSite(index));
        KERNEL::LatticeImpl::CalculateDensityMomentumFEq(
            hydroVars.f, hydroVars.density,
            hydroVars.momentum.x, hydroVars.momentum.y, hydroVars.momentum.z,
            hydroVars.velocity.x, hydroVars.velocity.y, hydroVars.velocity.z,
            hydroVars.f_eq.f
        );
        return hydroVars.velocity;
      }
      protected:
        //! Lattice data that holds the grid of population and forces
        LatticeData const & latticeData;
    };
  }

  template<class KERNEL>
    PhysicalVelocity interpolateVelocity(
        LatticeData const &_latdat, PhysicalPosition const &_center
        stencil::type _stencil = stencil::FOUR_POINT) {
      return interpolate(
          details::VelocityFromLatticeData<KERNEL>(_latdat),
          interpolation_iterator(_center, _stencil)
      );
    }
}} // hemelb::redblood
#endif
