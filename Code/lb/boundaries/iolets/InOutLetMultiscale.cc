// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "lb/boundaries/iolets/InOutLetMultiscale.h"
#include "configuration/SimConfig.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {
        void InOutLetMultiscale::DoIO(TiXmlElement *parent, bool isLoading, configuration::SimConfig* simConfig)
        {
          simConfig->DoIOForMultiscaleInOutlet(parent, isLoading, this);
        }
      }
    }
  }
}
