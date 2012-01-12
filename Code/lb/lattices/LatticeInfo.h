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
                      const Direction* inverseVectorIndicesIn);

          unsigned GetNumVectors() const;

          const util::Vector3D<int>& GetVector(unsigned index) const;

          unsigned GetInverseIndex(unsigned index) const;

        private:
          const unsigned numVectors;
          std::vector<util::Vector3D<int> > vectorSet;
          std::vector<Direction> inverseVectorIndices;
      };
    }
  }
}
#endif
