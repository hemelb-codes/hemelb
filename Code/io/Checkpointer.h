// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_CHECKPOINTER_H
#define HEMELB_IO_CHECKPOINTER_H

#include <filesystem>
#include <memory>

#include "io/PathManager.h"
#include "io/TimePattern.h"
#include "lb/SimulationState.h"
#include "net/IteratedAction.h"
#include "reporting/timers_fwd.h"

namespace hemelb::extraction { class LocalPropertyOutput; }

namespace hemelb::io {

    // Control the writing of a complete representation of a simulation state in order to restart.
    class Checkpointer : public net::IteratedAction {
        LatticeTimeStep period;

        // Path to each individual checkpoint's pattern, with a single "%d" for the time
        io::TimePattern out_dir_pattern;

        std::shared_ptr<lb::SimulationState const> simulationState;
        reporting::Timers& timers;

        // Delegate writing of the distribution array to property extraction
        std::unique_ptr<extraction::LocalPropertyOutput> distribution_writer;

    public:
        inline Checkpointer(
                LatticeTimeStep per,
                std::string_view pattern,
                std::shared_ptr<lb::SimulationState const> ss,
                reporting::Timers& t,
                std::unique_ptr<extraction::LocalPropertyOutput> dist_writer
        ) :
                period(per),
                out_dir_pattern(pattern),
                simulationState(std::move(ss)),
                timers(t),
                distribution_writer(std::move(dist_writer))
        {}

        ~Checkpointer() override;

        // Called every time step to give the option of writing.
        void EndIteration() override;

    private:
        [[nodiscard]] bool ShouldWrite() const;
        void Write() const;
    };
}

#endif