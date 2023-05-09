// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/FieldData.h"

#include "geometry/NeighbouringProcessor.h"
#include "geometry/neighbouring/NeighbouringDomain.h"
#include "lb/lattices/LatticeInfo.h"
#include "net/net.h"

namespace hemelb::geometry {
    FieldData::FieldData(std::shared_ptr <domain_type> d) :
            m_domain{d},
            m_currentDistributions(CalcDistSize(*d)),
            m_nextDistributions(CalcDistSize(*d)),
            m_force(d->GetLocalFluidSiteCount()),
            m_neighbouringFields{std::make_unique<neighbouring::NeighbouringFieldData>(d->neighbouringData)} {

    }

    std::size_t FieldData::CalcDistSize(Domain const &d) {
        return d.GetLocalFluidSiteCount() * d.latticeInfo.GetNumVectors() + 1 + d.totalSharedFs;
    }

    void FieldData::SendAndReceive(net::Net *net) {
        for (auto const &proc: GetDomain().neighbouringProcs) {
            // Request the receive into the appropriate bit of FOld.
            net->RequestReceive<distribn_t>(GetFOld(proc.FirstSharedDistribution),
                                            (int) proc.SharedDistributionCount,
                                            proc.Rank);
            // Request the send from the right bit of FNew.
            net->RequestSend<distribn_t>(GetFNew(proc.FirstSharedDistribution),
                                         (int) proc.SharedDistributionCount,
                                         proc.Rank);

        }
    }

    void FieldData::CopyReceived() {
        auto const &dom = GetDomain();
        // Copy the distribution functions received from the neighbouring
        // processors into the destination buffer "f_new".
        for (site_t i = 0; i < dom.totalSharedFs; i++) {
            *GetFNew(dom.streamingIndicesForReceivedDistributions[i]) =
                    *GetFOld(dom.neighbouringProcs[0].FirstSharedDistribution + i);
        }
    }

}