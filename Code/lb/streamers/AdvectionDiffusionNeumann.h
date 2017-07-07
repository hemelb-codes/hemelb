
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONNEUMANN_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONNEUMANN_H

#include "lb/streamers/AdvectionDiffusionStreamerTypeFactory.h"
#include "lb/streamers/AdvectionDiffusionNeumannDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"
#include "lb/streamers/GuoZhengShiDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<class CollisionType>
      struct AdvectionDiffusionNeumann
      {
          typedef StentWallStreamerTypeFactory<CollisionType, AdvectionDiffusionNeumannDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionNeumannSBB
      {
          typedef StentWallVesselWallStreamerTypeFactory<CollisionType, AdvectionDiffusionNeumannDelegate<CollisionType> ,
              SimpleBounceBackDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionNeumannBFL
      {
          typedef StentWallVesselWallStreamerTypeFactory<CollisionType, AdvectionDiffusionNeumannDelegate<CollisionType> ,
              BouzidiFirdaousLallemandDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionNeumannGZS
      {
          typedef StentWallVesselWallStreamerTypeFactory<CollisionType, AdvectionDiffusionNeumannDelegate<CollisionType> ,
              GuoZhengShiDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionNeumannIoletSBB
      {
          typedef StentWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionNeumannDelegate<CollisionType> ,
              SimpleBounceBackDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionNeumannIoletBFL
      {
          typedef StentWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionNeumannDelegate<CollisionType> ,
              BouzidiFirdaousLallemandDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionNeumannIoletGZS
      {
          typedef StentWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionNeumannDelegate<CollisionType> ,
              GuoZhengShiDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionNeumannSBBIoletWall
      {
          typedef StentWallVesselWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionNeumannDelegate<CollisionType> ,
              SimpleBounceBackDelegate<CollisionType> , SimpleBounceBackDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionNeumannBFLIoletWall
      {
          typedef StentWallVesselWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionNeumannDelegate<CollisionType> ,
              BouzidiFirdaousLallemandDelegate<CollisionType> , BouzidiFirdaousLallemandDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct AdvectionDiffusionNeumannGZSIoletWall
      {
          typedef StentWallVesselWallIoletStreamerTypeFactory<CollisionType, AdvectionDiffusionNeumannDelegate<CollisionType> ,
              GuoZhengShiDelegate<CollisionType> , GuoZhengShiDelegate<CollisionType> > Type;
      };


    }
  }
}

#endif /* HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONDIRICHLET_H */
