// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_DECOMPOSITION_BASICDECOMPOSITION_H
#define HEMELB_GEOMETRY_DECOMPOSITION_BASICDECOMPOSITION_H

#include <vector>

#include "units.h"
#include "util/Vector3D.h"

namespace hemelb {
    namespace geometry { class GmyReadResult; }
    namespace net { class MpiCommunicator; }

    namespace geometry::octree { class LookupTree; }
    namespace geometry::decomposition {
        class BasicDecomposition {
        public:
            using BlockLocation = util::Vector3D<site_t>;

            /**
             * Constructor to populate all fields necessary for a decomposition
             *
             * NOTE: We need the geometry (and its associated octree) in order to try to keep
             * contiguous blocks together, and to skip blocks with no fluid sites.
             *
             * @param geometry
             * @param communicator
             * @param fluidSitesOnEachBlock
             */
            BasicDecomposition(const GmyReadResult &geometry,
                               const net::MpiCommunicator &communicator);

            /**
             * Does a basic decomposition of the geometry without requiring any communication;
             * produces a vector of the processor assigned to each block.
             *
             * To make this fast just assign in octree order.
             *
             * @param rankForEachBlock A vector with the processor rank each block has been assigned to.
             *
             * @return a vector of the rank assigned to each non-solid block (i.e. leaf nodes on the tree).
             */
            std::vector<int>
            Decompose(octree::LookupTree const &blockTree, std::vector<proc_t> &procAssignedToEachBlock);

            /**
             * Validates that all cores have the same beliefs about which proc is to be assigned
             * to each proc by this decomposition.
             *
             * @param procAssignedToEachBlock This core's decomposition result.
             */
            void Validate(std::vector<proc_t> &procAssignedToEachBlock);

        private:

            const GmyReadResult &geometry; //! The geometry being decomposed.
            const net::MpiCommunicator &communicator; //! The communicator object being decomposed over.
        };
    }
}
#endif
