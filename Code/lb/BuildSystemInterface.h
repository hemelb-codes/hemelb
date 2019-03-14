
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_BUILDSYSTEMINTERFACE_H
#define HEMELB_LB_BUILDSYSTEMINTERFACE_H

#include "lb/kernels/Kernels.h"
#include "lb/kernels/momentBasis/MomentBases.h"
#include "lb/kernels/rheologyModels/RheologyModels.h"
#include "lb/collisions/Collisions.h"
#include "lb/streamers/Streamers.h"

namespace hemelb
{
  namespace lb
  {
    /**
     * The names of the following classes must correspond to options given for the CMake
     * HEMELB_KERNEL parameter.
     */
    template<class Lattice>
    class LBGK
    {
      public:
        typedef kernels::LBGK<Lattice> Type;
    };

    template<class Lattice>
    class ADVECTIONDIFFUSIONLBGK
    {
      public:
        typedef kernels::AdvectionDiffusionLBGK<Lattice> Type;
    };

    /**
     * The entropic implementation by Ansumali et al.
     */
    template<class Lattice>
    class EntropicAnsumali
    {
      public:
        typedef kernels::EntropicAnsumali<Lattice> Type;
    };

    /**
     * The entropic implementation by Chikatamarla et al.
     */
    template<class Lattice>
    class EntropicChik
    {
      public:
        typedef kernels::EntropicChik<Lattice> Type;
    };

    /**
     * MRT currently we only have DHumieres implementation, on D3Q15 and D3Q19 lattices.
     */
    template<class Lattice>
    class MRT
    {
    };

    template<class Lattice>
    class ADVECTIONDIFFUSIONMRT
    {
    };

    template<>
    class MRT<lattices::D3Q15>
    {
      public:
        typedef kernels::MRT<kernels::momentBasis::DHumieresD3Q15MRTBasis> Type;
    };

    template<>
    class ADVECTIONDIFFUSIONMRT<lattices::D3Q15>
    {
      public:
        typedef kernels::AdvectionDiffusionMRT<kernels::momentBasis::DHumieresD3Q15MRTBasis> Type;
    };

    template<>
    class MRT<lattices::D3Q19>
    {
      public:
        typedef kernels::MRT<kernels::momentBasis::DHumieresD3Q19MRTBasis> Type;
    };

    template<>
    class ADVECTIONDIFFUSIONMRT<lattices::D3Q19>
    {
      public:
        typedef kernels::AdvectionDiffusionMRT<kernels::momentBasis::DHumieresD3Q19MRTBasis> Type;
    };

    template<class Lattice>
    class TRT
    {
      public:
        typedef kernels::TRT<Lattice> Type;
    };

    /**
     * Non-Newtonian kernel with Carreau-Yasuda rheology model.
     */
    template<class Lattice>
    class NNCY
    {
      public:
        typedef kernels::LBGKNN<kernels::rheologyModels::CarreauYasudaRheologyModelHumanFit, Lattice> Type;
    };

    /**
     * Non-Newtonian kernel with Carreau-Yasuda rheology model fitted to experimental data on murine blood viscosity.
     */
    template<class Lattice>
    class NNCYMOUSE
    {
      public:
        typedef kernels::LBGKNN<kernels::rheologyModels::CarreauYasudaRheologyModelMouseFit, Lattice> Type;
    };

    /**
     * Non-Newtonian kernel with Casson rheology model.
     */
    template<class Lattice>
    class NNC
    {
      public:
        typedef kernels::LBGKNN<kernels::rheologyModels::CassonRheologyModel, Lattice> Type;
    };

    /**
     * Non-Newtonian kernel with truncated power law rheology model.
     */
    template<class Lattice>
    class NNTPL
    {
      public:
        typedef kernels::LBGKNN<kernels::rheologyModels::TruncatedPowerLawRheologyModel, Lattice> Type;
    };

    /**
     * The following classes have names corresponding to the options given in the build system for
     * HEMELB_WALL_BOUNDARY
     */
    /**
     * The Bouzidi-Firdaous-Lallemand interpolation-based boundary condition.
     */
    template<class Collision>
    class BFL
    {
      public:
        typedef typename streamers::BouzidiFirdaousLallemand<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLCoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLIOLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLIolet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLBFLIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallBFLIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLBFLIOLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallBFLIoletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLBFLIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallBFLIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLBFLIOLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallBFLIoletCoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLSBBIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallSBBIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLSBBIOLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallSBBIoletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLSBBIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallSBBIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLSBBIOLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallSBBIoletCoatingFlux<Collision>::Type Type;
    };
    /**
     * The Guo Zheng and Shi mode-extrapolation boundary condition.
     */
    template<class Collision>
    class GZS
    {
      public:
        typedef typename streamers::GuoZhengShi<Collision>::Type Type;
    };
    /**
     * The simple bounce back boundary condition.
     */
    template<class Collision>
    class SIMPLEBOUNCEBACK
    {
      public:
        typedef typename streamers::SimpleBounceBack<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBCoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBIOLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBIolet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLBFLIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallBFLIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLBFLIOLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallBFLIoletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLBFLIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallBFLIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLBFLIOLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallBFLIoletCoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLSBBIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallSBBIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLSBBIOLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallSBBIoletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLSBBIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallSBBIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLSBBIOLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallSBBIoletCoatingFlux<Collision>::Type Type;
    };
    /**
     * The Junk & Yang 2005 boundary condition.
     */
    template<class Collision>
    class JUNKYANG
    {
      public:
        typedef typename streamers::JunkYang<Collision>::Type Type;
    };

    /**
     * The following classes have names corresponding to the options given in the build system for
     * HEMELB_INLET_BOUNDARY / HEMELB_OUTLET_BOUNDARY
     */

    /**
     * Our zeroth-order phantom site BC for iolets
     */
    template<class Collision>
    class NASHZEROTHORDERPRESSUREIOLET
    {
      public:
        typedef typename streamers::NashZerothOrderPressureIolet<Collision>::Type Type;
    };
    /**
     * The inlet/outlet condition based on Ladd's modified bounce-back on
     * links.
     */
    template<class Collision>
    struct LADDIOLET
    {
        typedef typename streamers::LaddIolet<Collision>::Type Type;
    };

    /**
     * The following classes have names corresponding to the options given in the build system for
     * HEMELB_WALL_INLET_BOUNDARY / HEMELB_WALL_OUTLET_BOUNDARY
     */
    /**
     * Nash in/outlet + SBB
     */
    template<class Collision>
    class NASHZEROTHORDERPRESSURESBB
    {
      public:
        typedef typename streamers::NashZerothOrderPressureIoletSBB<Collision>::Type Type;
    };

    /**
     * Ladd in/outlet + SBB
     */
    template<class Collision>
    struct LADDIOLETSBB
    {
        typedef typename streamers::LaddIoletSBB<Collision>::Type Type;
    };

    /**
     * Nash in/outlet + BFL
     */
    template<class Collision>
    class NASHZEROTHORDERPRESSUREBFL
    {
      public:
        typedef typename streamers::NashZerothOrderPressureIoletBFL<Collision>::Type Type;
    };

    /**
     * Ladd in/outlet + BFL
     */
    template<class Collision>
    struct LADDIOLETBFL
    {
        typedef typename streamers::LaddIoletBFL<Collision>::Type Type;
    };
    /**
     * Nash in/outlet + GZS
     */
    template<class Collision>
    class NASHZEROTHORDERPRESSUREGZS
    {
      public:
        typedef typename streamers::NashZerothOrderPressureIoletGZS<Collision>::Type Type;
    };

    /**
     * Ladd in/outlet + GZS
     */
    template<class Collision>
    struct LADDIOLETGZS
    {
        typedef typename streamers::LaddIoletGZS<Collision>::Type Type;
    };

    /**
     * Nash/Krueger in/outlet + SBB
     */
    template<class Collision>
    struct VIRTUALSITEIOLETSBB
    {
        typedef typename streamers::VirtualSiteIolet<Collision> Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONOUTFLOWIOLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionOutflowIolet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWIOLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowIoletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWIOLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowIoletCoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWIOLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowIoletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWIOLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowIoletCoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONINLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionInletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLINLETDIRICHLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallInletDirichletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLINLETDIRICHLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallInletDirichletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLINLETDIRICHLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallInletDirichletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLINLETDIRICHLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallInletDirichletCoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLINLETDIRICHLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallInletDirichletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLINLETDIRICHLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallInletDirichletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLINLETDIRICHLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallInletDirichletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLINLETDIRICHLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallInletDirichletCoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONOUTFLOWBOUNCEBACKIOLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionOutflowBounceBackIolet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWBOUNCEBACKIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowBounceBackIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWBOUNCEBACKIOLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowBounceBackIoletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWBOUNCEBACKIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowBounceBackIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONSBBWALLOUTFLOWBOUNCEBACKIOLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionSBBWallOutflowBounceBackIoletCoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWBOUNCEBACKIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowBounceBackIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWBOUNCEBACKIOLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowBounceBackIoletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWBOUNCEBACKIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowBounceBackIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONBFLWALLOUTFLOWBOUNCEBACKIOLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionBFLWallOutflowBounceBackIoletCoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWADIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWADirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWACOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWACoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWANEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWANeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWACOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWACoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWAWALLINLETDIRICHLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWAWallInletDirichletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWAWALLINLETDIRICHLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWAWallInletDirichletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWAWALLINLETDIRICHLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWAWallInletDirichletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWAWALLINLETDIRICHLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWAWallInletDirichletCoatingFlux<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWAWALLOUTFLOWBOUNCEBACKIOLETDIRICHLET
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWAWallOutflowBounceBackIoletDirichlet<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWAWALLOUTFLOWBOUNCEBACKIOLETCOATINGCONCENTRATION
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWAWallOutflowBounceBackIoletCoatingConcentration<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWAWALLOUTFLOWBOUNCEBACKIOLETNEUMANN
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWAWallOutflowBounceBackIoletNeumann<Collision>::Type Type;
    };

    template<class Collision>
    class ADVECTIONDIFFUSIONVWAWALLOUTFLOWBOUNCEBACKIOLETCOATINGFLUX
    {
      public:
        typedef typename streamers::AdvectionDiffusionVWAWallOutflowBounceBackIoletCoatingFlux<Collision>::Type Type;
    };
  }
}

#endif /* HEMELB_LB_BUILDSYSTEMINTERFACE_H */
