// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_VELOCITYINTERPOLATION_H
#define HEMELB_REDBLOOD_VELOCITYINTERPOLATION_H

#include <iostream>

#include "units.h"
#include "Exception.h"
#include "redblood/Interpolation.h"
#include "redblood/stencil.h"
#include "lb/kernels/GuoForcingLBGK.h"

#include <vector>
#include <boost/shared_array.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/bool.hpp>

namespace hemelb
{
  namespace redblood
  {
    //! \brief Returns the interpolated velocity at given point
    //! \param[in] center: off lattice position for which to interpolate.
    //! \param[in] latDat: lattice data
    //! \param[in] stencil: stencil to use for interpolation
    //! \tparam KERNEL: Needed to compute velocity
    template<class KERNEL>
    LatticeVelocity interpolateVelocity(geometry::LatticeData const &latDat,
                                        LatticePosition const &center);

    namespace details
    {
      // Type traits to specialize for kernels with forces
      template<class KERNEL>
      struct HasForce : public boost::mpl::false_
      {
      };
      template<class LATTICE>
      struct HasForce<lb::kernels::GuoForcingLBGK<LATTICE> > : public boost::mpl::true_
      {
      };

      // Computes velocity for a given index on the lattice
      template<class KERNEL>
      struct VelocityFromLatticeData
      {
          VelocityFromLatticeData(geometry::LatticeData const &latDat) :
              latticeData(latDat)
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

      template<class KERNEL>
      LatticeVelocity VelocityFromLatticeData<KERNEL>::computeVelocity(
          site_t index, boost::mpl::true_ const &) const
      {
        typedef typename KERNEL::LatticeType LatticeType;
        LatticeVelocity result;
        distribn_t density;
        geometry::Site<geometry::LatticeData const> site(latticeData.GetSite(index));
        LatticeForceVector const &force(site.GetForce());
#ifdef HEMELB_USE_KRUEGER_ORDERING
        // Use distribution functions at the beginning of the previous timestep (stored in
        // FNew after the swap at the end of the timestep) in the IBM velocity interpolation.
        // Follows approach in Timm's code
        auto const fDistribution = latticeData.GetFNew(index * LatticeType::NUMVECTORS);
#else
        auto const fDistribution = site.GetFOld<LatticeType>();
#endif
        LatticeType::CalculateDensityAndMomentum(fDistribution,
                                                 force[0],
                                                 force[1],
                                                 force[2],
                                                 density,
                                                 result.x,
                                                 result.y,
                                                 result.z);
        return result / density;
      }

      template<class KERNEL>
      LatticeVelocity VelocityFromLatticeData<KERNEL>::computeVelocity(
          site_t index, boost::mpl::false_ const &) const
      {
        typedef typename KERNEL::LatticeType LatticeType;
        LatticeVelocity result(0, 0, 0);
        distribn_t density;
#ifdef HEMELB_USE_KRUEGER_ORDERING
        // Use distribution functions at the beginning of the previous timestep (stored in
        // FNew after the swap at the end of the timestep) in the IBM velocity interpolation.
        // Follows approach in Timm's code
        auto const fDistribution = latticeData.GetFNew(index * LatticeType::NUMVECTORS);
#else
        auto const fDistribution = latticeData.GetSite(index).template GetFOld<LatticeType>();
#endif
        LatticeType::CalculateDensityAndMomentum(fDistribution,
                                                 density,
                                                 result.x,
                                                 result.y,
                                                 result.z);
        return result / density;
      }
    }

    template<class KERNEL, class STENCIL>
    LatticeVelocity interpolateVelocity(geometry::LatticeData const &latDat,
                                        LatticePosition const &center)
    {
      auto iterator = interpolationIterator<STENCIL>(center);
      // Computes velocity for a given site index
      // Branches to one or another function depending on whether forces are available (since
      // velocity depends on forces)
      auto const gridfunc = details::VelocityFromLatticeData<KERNEL>(latDat);
      LatticeForceVector result(0, 0, 0);
      for (; iterator.IsValid(); ++iterator)
      {
        proc_t procid;
        site_t siteid;
        if (latDat.GetContiguousSiteId(*iterator, procid, siteid))
        {
          result += gridfunc(siteid) * iterator.weight();
        }
      }
      return result;
    }
  }
} // hemelb::redblood
#endif
