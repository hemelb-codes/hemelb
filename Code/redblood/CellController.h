//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELL_CONTROLLER_H
#define HEMELB_UNITTESTS_REDBLOOD_CELL_CONTROLLER_H

#include "net/net.h"
#include "net/IteratedAction.h"
#include "redblood/CellArmy.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace redblood
  {
    //! \brief Federates the cells together so we can apply ops simultaneously
    template<class KERNEL>
      class CellController : public CellArmy<KERNEL>, public net::IteratedAction
    {
      public:
        using CellArmy<KERNEL>::CellArmy;

        void RequestComms() override
        {
          using namespace log;
          Logger::Log<Info, Singleton>("Cell insertion");
          CellArmy<KERNEL>::CallCellInsertion();
          Logger::Log<Info, Singleton>("Fluid interaction with cells");
          CellArmy<KERNEL>::Fluid2CellInteractions();
        }
        void EndIteration() override
        {
          using namespace log;
          Logger::Log<Info, Singleton>("Cell interaction with fluid");
          CellArmy<KERNEL>::Cell2FluidInteractions(); 
          Logger::Log<Info, Singleton>("Removed cells that have reached outlets");
          CellArmy<KERNEL>::CellRemoval(); 
        }
    };
  }
}

#endif
