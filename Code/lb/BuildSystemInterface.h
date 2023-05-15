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

namespace hemelb::lb
{
    /**
     * The names of the following classes must correspond to options given for the CMake
     * HEMELB_KERNEL parameter.
     */
    template<class Lattice>
    class LBGK
    {
      public:
        using Type = kernels::LBGK<Lattice>;
    };

    /**
     * LBGK with GuoForcing
     */
    template<class LATTICE>
    class GuoForcingLBGK
    {
      public:
        using Type = lb::kernels::GuoForcingLBGK<LATTICE>;
    };

    /**
     * The entropic implementation by Ansumali et al.
     */
    template<class Lattice>
    class EntropicAnsumali
    {
      public:
        using Type = kernels::EntropicAnsumali<Lattice>;
    };

    /**
     * The entropic implementation by Chikatamarla et al.
     */
    template<class Lattice>
    class EntropicChik
    {
      public:
        using Type = kernels::EntropicChik<Lattice>;
    };

    /**
     * MRT currently we only have DHumieres implementation, on D3Q15 and D3Q19 lattices.
     */
    template<class Lattice>
    class MRT
    {
    };

    template<>
    class MRT<D3Q15>
    {
      public:
        using Type = kernels::MRT<kernels::momentBasis::DHumieresD3Q15MRTBasis>;
    };

    template<>
    class MRT<D3Q19>
    {
      public:
        using Type = kernels::MRT<kernels::momentBasis::DHumieresD3Q19MRTBasis>;
    };

    template<class Lattice>
    class TRT
    {
      public:
        using Type = kernels::TRT<Lattice>;
    };

    /**
     * Non-Newtonian kernel with Carreau-Yasuda rheology model.
     */
    template<class Lattice>
    class NNCY
    {
      public:
        using Type = kernels::LBGKNN<kernels::rheologyModels::CarreauYasudaRheologyModelHumanFit, Lattice>;
    };

    /**
     * Non-Newtonian kernel with Carreau-Yasuda rheology model fitted to experimental data on murine blood viscosity.
     */
    template<class Lattice>
    class NNCYMOUSE
    {
      public:
        using Type = kernels::LBGKNN<kernels::rheologyModels::CarreauYasudaRheologyModelMouseFit, Lattice>;
    };

    /**
     * Non-Newtonian kernel with Casson rheology model.
     */
    template<class Lattice>
    class NNC
    {
      public:
        using Type = kernels::LBGKNN<kernels::rheologyModels::CassonRheologyModel, Lattice>;
    };

    /**
     * Non-Newtonian kernel with truncated power law rheology model.
     */
    template<class Lattice>
    class NNTPL
    {
      public:
        using Type = kernels::LBGKNN<kernels::rheologyModels::TruncatedPowerLawRheologyModel, Lattice>;
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
        using Type = typename streamers::BouzidiFirdaousLallemand<Collision>::Type;
    };
    /**
     * The Guo Zheng and Shi mode-extrapolation boundary condition.
     */
    template<class Collision>
    class GZS
    {
      public:
        using Type = typename streamers::GuoZhengShi<Collision>::Type;
    };
    /**
     * The simple bounce back boundary condition.
     */
    template<class Collision>
    class SIMPLEBOUNCEBACK
    {
      public:
        using Type = typename streamers::SimpleBounceBack<Collision>::Type;
    };
    /**
     * The Junk & Yang 2005 boundary condition.
     */
    template<class Collision>
    class JUNKYANG
    {
      public:
        using Type = typename streamers::JunkYang<Collision>::Type;
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
        using Type = typename streamers::NashZerothOrderPressureIolet<Collision>::Type;
    };
    /**
     * The inlet/outlet condition based on Ladd's modified bounce-back on
     * links.
     */
    template<class Collision>
    struct LADDIOLET
    {
        using Type = typename streamers::LaddIolet<Collision>::Type;
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
        using Type = typename streamers::NashZerothOrderPressureIoletSBB<Collision>::Type;
    };

    /**
     * Ladd in/outlet + SBB
     */
    template<class Collision>
    struct LADDIOLETSBB
    {
        using Type = typename streamers::LaddIoletSBB<Collision>::Type;
    };

    /**
     * Nash in/outlet + BFL
     */
    template<class Collision>
    class NASHZEROTHORDERPRESSUREBFL
    {
      public:
        using Type = typename streamers::NashZerothOrderPressureIoletBFL<Collision>::Type;
    };

    /**
     * Ladd in/outlet + BFL
     */
    template<class Collision>
    struct LADDIOLETBFL
    {
        using Type = typename streamers::LaddIoletBFL<Collision>::Type;
    };
    /**
     * Nash in/outlet + GZS
     */
    template<class Collision>
    class NASHZEROTHORDERPRESSUREGZS
    {
      public:
        using Type = typename streamers::NashZerothOrderPressureIoletGZS<Collision>::Type;
    };

    /**
     * Ladd in/outlet + GZS
     */
    template<class Collision>
    struct LADDIOLETGZS
    {
        using Type = typename streamers::LaddIoletGZS<Collision>::Type;
    };

    /**
     * Nash/Krueger in/outlet + SBB
     */
    template<class Collision>
    struct VIRTUALSITEIOLETSBB
    {
        using Type = typename streamers::VirtualSiteIolet<Collision>;
    };
  }

#endif /* HEMELB_LB_BUILDSYSTEMINTERFACE_H */
