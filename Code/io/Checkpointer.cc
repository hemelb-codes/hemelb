// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/Checkpointer.h"

#include "extraction/LocalPropertyOutput.h"
#include "reporting/Timers.h"

namespace hemelb::io {
    namespace fs = std::filesystem;

    Checkpointer::~Checkpointer() = default;

    void Checkpointer::EndIteration() {
        timers.writeCheckpoint().Start();
        if (ShouldWrite())
            Write();
        timers.writeCheckpoint().Stop();
    }

    bool Checkpointer::ShouldWrite() const {
        return (simulationState->GetTimeStep() % period) == 0;
    }

    void Checkpointer::Write() const {
        auto dir = out_dir_pattern.Format(simulationState->GetTimeStep(), simulationState->GetEndTimeStep());
        fs::create_directory(dir);

        distribution_writer->Write(simulationState->GetTimeStep(), simulationState->GetEndTimeStep());

    }
}