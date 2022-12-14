// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LATTICES_LATTICEINFO_H
#define HEMELB_LB_LATTICES_LATTICEINFO_H

#include <array>
#include "Exception.h"
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb::lb::lattices
{
    class LatticeInfo
    {
    public:
        constexpr LatticeInfo(
                std::size_t numberOfVectors,
                const util::Vector3D<int>* vectors,
                const Direction* inverseVectorIndicesIn
        ) : numVectors(numberOfVectors)
        {
            if (numberOfVectors > MAX_Q)
                throw Exception() << "LatticeInfo only supports velocity sets with a maximum of " << MAX_Q << " velocities.";

            std::copy(vectors, vectors + numberOfVectors, vectorSet.begin());
            std::copy(inverseVectorIndicesIn, inverseVectorIndicesIn + numberOfVectors, inverseVectorIndices.begin());
        }

        constexpr unsigned GetNumVectors() const
        {
            return numVectors;
        }

        constexpr const util::Vector3D<int>& GetVector(unsigned index) const
        {
            return vectorSet[index];
        }

        constexpr unsigned GetInverseIndex(unsigned index) const
        {
            return inverseVectorIndices[index];
        }

    private:
        // Max possible number to help with constexpr-ness
        static constexpr std::size_t MAX_Q = 27;
        // Actual number
        unsigned numVectors;
        // Make storage for maximum number
        std::array<util::Vector3D<int>, MAX_Q> vectorSet;
        std::array<Direction, MAX_Q> inverseVectorIndices;
    };
}
#endif
