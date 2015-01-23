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

#include <iostream>

#include "units.h"
#include "Exception.h"
#include "redblood/Interpolation.h"
#include "redblood/stencil.h"
#include "lb/kernels/GuoForcingLBGK.h"

#include <vector>
#include <cassert>
#include <boost/shared_array.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/bool.hpp>

namespace hemelb
{
  namespace redblood
  {
    //! \brief Returns the interpolated velocity at given point
    //! \param[in] center: off lattice position for which to interpolate.
    //! \param[in] latdat: lattice data
    //! \param[in] stencil: stencil to use for interpolation
    //! \tparam KERNEL: Needed to compute velocity
    template <class KERNEL>
    PhysicalVelocity interpolateVelocity(PhysicalPosition const &center,
                                         geometry::LatticeData const &latdat,
                                         stencil::types stencil = stencil::FOUR_POINT);

    namespace details
    {
      // Type traits to specialize for kernels with forces
      template <class KERNEL>
      struct HasForce : public boost::mpl::false_
      {
      };
      template <class LATTICE>
      struct HasForce<lb::kernels::GuoForcingLBGK<LATTICE> > : public boost::mpl::true_
      {
      };

      // Computes velocity for a given index on the lattice
      template <class KERNEL>
      struct VelocityFromLatticeData
      {
        VelocityFromLatticeData(geometry::LatticeData const &latdata) : latticeData(latdata)
        {
        }
        // Computes velocity for given KERNEL and LatticeData
        LatticeVelocity operator()(LatticeVector const &indices) const
        {
          return operator()(latticeData.GetContiguousSiteId(indices));
        }
        LatticeVelocity operator()(size_t index) const
        {
          return computeVelocity(index, HasForce<KERNEL>());
        }

        protected:
        //! Implementation with forces correction
        LatticeVelocity computeVelocity(site_t index, boost::mpl::true_ const &) const;
        //! Implementation without forces correction
        LatticeVelocity computeVelocity(site_t index, boost::mpl::false_ const &) const;
        //! Lattice data that holds the grid of population and forces
        geometry::LatticeData const &latticeData;
      };

      template <class KERNEL>
      LatticeVelocity VelocityFromLatticeData<KERNEL>::computeVelocity(
        site_t index, boost::mpl::true_ const &) const
      {
        typedef typename KERNEL::LatticeType LatticeType;
        LatticeVelocity result;
        distribn_t density;
        geometry::Site<geometry::LatticeData const> site(latticeData.GetSite(index));
        LatticeForceVector const &force(site.GetForce());
        LatticeType::CalculateDensityAndMomentum(site.GetFOld<LatticeType>(), force[0], force[1],
                                                 force[2], density, result.x, result.y, result.z);
        return result / density;
      }

      template <class KERNEL>
      LatticeVelocity VelocityFromLatticeData<KERNEL>::computeVelocity(
        site_t index, boost::mpl::false_ const &) const
      {
        typedef typename KERNEL::LatticeType LatticeType;
        LatticeVelocity result(0, 0, 0);
        distribn_t density;
        LatticeType::CalculateDensityAndMomentum(
          latticeData.GetSite(index).template GetFOld<LatticeType>(), density, result.x, result.y,
          result.z);
        return result / density;
      }
    }

    template <class KERNEL>
    PhysicalVelocity interpolateVelocity(geometry::LatticeData const &latdat,
                                         PhysicalPosition const &center,
                                         stencil::types stencil = stencil::FOUR_POINT)
    {
      return interpolate(details::VelocityFromLatticeData<KERNEL>(latdat),
                         interpolationIterator(center, stencil));
    }
  }
}  // hemelb::redblood
#endif
