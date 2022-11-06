// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_PROPERTYOUTPUTFILE_H
#define HEMELB_EXTRACTION_PROPERTYOUTPUTFILE_H

#include <filesystem>
#include <variant>
#include <vector>

#include "util/clone_ptr.h"
#include "extraction/GeometrySelector.h"
#include "extraction/OutputField.h"

namespace hemelb::extraction
{
  // Tag types for the file timestep mode.
  struct multi_timestep_file {};
  struct single_timestep_files {};

  // Multiple timesteps per file first so it will be default (as this
  // is the old behaviour).
  using file_timestep_mode = std::variant<multi_timestep_file, single_timestep_files>;

  struct PropertyOutputFile
  {
    std::filesystem::path filename;
    unsigned long frequency;
    util::clone_ptr<GeometrySelector> geometry;
    std::vector<OutputField> fields;
    file_timestep_mode ts_mode;
  };
}

#endif // HEMELB_EXTRACTION_PROPERTYOUTPUTFILE_H
