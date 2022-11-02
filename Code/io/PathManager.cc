// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/PathManager.h"

#include "Exception.h"
#include "configuration/CommandLine.h"

namespace fs = std::filesystem;
namespace hemelb::io
{

    PathManager::PathManager(const configuration::CommandLine & commandLine, const bool & io,
                             const int & processorCount) :
            options(commandLine), doIo(io)
    {
        inputFile = fs::absolute(options.GetInputFile());
        outputDir = fs::absolute(options.GetOutputDir());

        extractionDir = outputDir / "Extracted";
        colloidFile = outputDir / "ColloidOutput.xdr";
        rbcDir = outputDir / "Cells";

        if (doIo)
        {
            if (fs::exists(outputDir))
                throw Exception() << "Output directory '" << outputDir << "' already exists.";

            fs::create_directory(outputDir);
            fs::create_directory(extractionDir);
            fs::create_directory(rbcDir);
        }
    }

    const fs::path& PathManager::GetInputFile() const
    {
      return inputFile;
    }
    const fs::path& PathManager::GetColloidPath() const
    {
      return colloidFile;
    }
    const fs::path& PathManager::GetReportPath() const
    {
      return outputDir;
    }

    const fs::path& PathManager::GetDataExtractionPath() const
    {
      return extractionDir;
    }

    void PathManager::SaveConfiguration(configuration::SimConfig * const simConfig) const
    {
      if (doIo)
      {
        //simConfig->Save(outputDir + "/" + configLeafName);
      }
    }

    fs::path PathManager::GetRBCOutputPathWithSubdir(std::string const& subdirectoryName) const
    {
      auto rbcSubdir = rbcDir / subdirectoryName;

      if (doIo)
      {
        if (fs::exists(rbcSubdir))
          throw Exception() << "Output directory '" << rbcSubdir << "' already exists.";

        std::error_code ec;
        if (!fs::create_directory(rbcSubdir, rbcDir, ec))
        if (ec)
        {
          throw Exception() << "Output directory '" << rbcSubdir << "' could not be created because: "
                            << ec.message();
        }
      }

      return rbcSubdir;
    }

}

