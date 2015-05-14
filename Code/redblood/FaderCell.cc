//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/FaderCell.h"
#include "redblood/FlowExtension.h"

namespace hemelb
{
  namespace redblood
  {
    PhysicalEnergy FaderCell::operator()() const
    {
      auto const barycenter = wrappee->GetBarycenter();
      auto const energy = wrappee->Energy();
      for (auto const& extension : *iolets)
      {
        auto const weight = linearWeight(extension, barycenter);
        if (weight > 1e-12)
        {
          return energy * weight;
        }
      }
      return energy;
    }

    PhysicalEnergy FaderCell::operator()(std::vector<LatticeForceVector> &forces) const
    {
      auto const barycenter = wrappee->GetBarycenter();
      auto const energy = wrappee->Energy(forces);
      for (auto const& extension : *iolets)
      {
        auto const weight = linearWeight(extension, barycenter);
        if (weight > 1e-12)
        {
          for (auto& force : forces)
          {
            force *= weight;
          }
          return energy * weight;
        }
      }
      return energy;
    }
  }
} // hemelb::redblood
