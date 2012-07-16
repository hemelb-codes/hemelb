// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_LATTICES_LATTICEINFO_H
#define HEMELB_LB_LATTICES_LATTICEINFO_H

#include <vector>
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      class LatticeInfo
      {
        public:
          inline LatticeInfo(unsigned numberOfVectors,
                             const util::Vector3D<int>* vectors,
                             const Direction* inverseVectorIndicesIn) :
              numVectors(numberOfVectors), vectorSet(), inverseVectorIndices()
          {
            for (Direction direction = 0; direction < numberOfVectors; ++direction)
            {
              vectorSet.push_back(util::Vector3D<int>(vectors[direction]));
              inverseVectorIndices.push_back(inverseVectorIndicesIn[direction]);
            }
          }

          inline unsigned GetNumVectors() const
          {
            return numVectors;
          }

          inline const util::Vector3D<int>& GetVector(unsigned index) const
          {
            return vectorSet[index];
          }

          inline unsigned GetInverseIndex(unsigned index) const
          {
            return inverseVectorIndices[index];
          }

        private:
          const unsigned numVectors;
          std::vector<util::Vector3D<int> > vectorSet;
          std::vector<Direction> inverseVectorIndices;
      };
    }
  }
}
#endif
