// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "lb/lattices/D3Q27.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      template<>
      LatticeInfo* Lattice<D3Q27>::singletonInfo = nullptr;

      const int D3Q27::CX[] =
          { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
      const int D3Q27::CY[] =
          { 0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1 };
      const int D3Q27::CZ[] =
          { 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1 };
      
      const distribn_t  D3Q27::CXD[] = 
          { 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0,  -1.0, 1.0, -1.0, 0.0, 0.0, 
            0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 };
      const distribn_t  D3Q27::CYD[] = 
          { 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0,
           -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0 };
      const distribn_t  D3Q27::CZD[] = 
         { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,  -1.0, 1.0, 1.0, -1.0, -1.0,
          1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0 };
      
      const int* D3Q27::discreteVelocityVectors[] = { CX, CY, CZ };

      const distribn_t D3Q27::EQMWEIGHTS[] = { 8.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0,
                                               2.0 / 27.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
                                               1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
                                               1.0 / 54.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0
                                                   / 216.0,
                                               1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0 };

      const Direction D3Q27::INVERSEDIRECTIONS[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17,
                                                     20, 19, 22, 21, 24, 23, 26, 25 };
    }
  }
}
