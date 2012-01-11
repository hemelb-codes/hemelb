#include "lb/lattices/LatticeInfo.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      LatticeInfo::LatticeInfo(unsigned numberOfVectors,
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

      unsigned LatticeInfo::GetNumVectors() const
      {
        return numVectors;
      }

      const util::Vector3D<int>& LatticeInfo::GetVector(unsigned index) const
      {
        return vectorSet[index];
      }

      unsigned LatticeInfo::GetInverseIndex(unsigned index) const
      {
        return inverseVectorIndices[index];
      }

    }
  }
}
