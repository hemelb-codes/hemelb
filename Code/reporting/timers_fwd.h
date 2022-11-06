// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REPORTING_TIMERS_FWD_H
#define HEMELB_REPORTING_TIMERS_FWD_H

namespace hemelb::reporting {
    template<class>
    class TimerBase;
    template<class ClockPolicy, class CommsPolicy>
    class TimersBase;

    class MPICommsPolicy;
    class HemeLBClockPolicy;

    using Timer = TimerBase<HemeLBClockPolicy>;
    using Timers = TimersBase<HemeLBClockPolicy, MPICommsPolicy>;

}
#endif