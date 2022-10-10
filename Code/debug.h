// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_DEBUG_H
#define HEMELB_DEBUG_H

namespace hemelb::net { class MpiCommunicator; }

namespace hemelb::debug {
    class Debugger;

    // Convenient interface to custom small-scale parallel debug.
    void Init(bool active, const char * executable, const net::MpiCommunicator& comm);
    void Break();
}

#endif