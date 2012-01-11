#ifndef HEMELB_LB_LATTICES_LATTICEINFO_H
#define HEMELB_LB_LATTICES_LATTICEINFO_H

#include "util/Vector3D.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      class LatticeInfo
      {
        public:
          LatticeInfo(unsigned numberOfVectors,
                      const util::Vector3D<int>* vectors,
                      const Direction* inverseVectorIndices) :
            NumVectors(numberOfVectors)
          {
            VectorSet = new util::Vector3D<int>[numberOfVectors];
            InverseVectorIndices = new unsigned[numberOfVectors];

            for (Direction direction = 0; direction < numberOfVectors; ++direction)
            {
              VectorSet[direction] = util::Vector3D<int>(vectors[direction]);
              InverseVectorIndices[direction] = inverseVectorIndices[direction];
            }
          }

          ~LatticeInfo()
          {
            delete[] VectorSet;
            delete[] InverseVectorIndices;
          }

          unsigned GetNumVectors() const
          {
            return NumVectors;
          }

          const util::Vector3D<int>& GetVector(unsigned index) const
          {
            return VectorSet[index];
          }

          const unsigned GetInverseIndex(unsigned index) const
          {
            return InverseVectorIndices[index];
          }

        private:
          const unsigned NumVectors;
          util::Vector3D<int>* VectorSet;
          Direction* InverseVectorIndices;
      };
    }
  }
}
#endif
