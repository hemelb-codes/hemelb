// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_LINKSTREAMER_H
#define HEMELB_LB_STREAMERS_LINKSTREAMER_H

#include "geometry/FieldData.h"

namespace hemelb::lb { class LbmParameters; }
namespace hemelb::lb::streamers
{
      /**
       * Concept for Streamer delegates.
       *
       * Sets out the interface that streamers should implement and typedefs
       * that they should make.
       */
      template <typename T>
      concept LinkStreamer = requires(
              T& linkStreamer,
              LbmParameters const* lbmParams,
              geometry::FieldData& data,
              geometry::Site<geometry::FieldData> const& site,
              HydroVars<typename T::CollisionType::CKernel>& hydroVars,
              Direction d
      ) {
          typename T::CollisionType;
          typename T::LatticeType;
          { linkStreamer.StreamLink(lbmParams, data, site, hydroVars, d) };
          { linkStreamer.PostStepLink(data, site, d) };
      };
}

#endif
