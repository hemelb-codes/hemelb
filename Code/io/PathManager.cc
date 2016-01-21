
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/PathManager.h"
#include "Exception.h"

namespace hemelb
{
  namespace io
  {
    PathManager::PathManager(const configuration::CommandLine & commandLine,
                             const bool & io,
                             const int & processorCount) :
      options(commandLine), doIo(io)
    {

      inputFile = options.GetInputFile();
      outputDir = options.GetOutputDir();

      GuessOutputDir();

      imageDirectory = outputDir + "/Images/";
      dataPath = outputDir + "/Extracted/";
      colloidFile = outputDir + "/ColloidOutput.xdr";

      if (doIo)
      {
        if (hemelb::util::DoesDirectoryExist(outputDir.c_str()))
          throw Exception() << "Output directory '"<< outputDir <<"' already exists.";

        hemelb::util::MakeDirAllRXW(outputDir);
        hemelb::util::MakeDirAllRXW(imageDirectory);
        hemelb::util::MakeDirAllRXW(dataPath);
        reportName = outputDir;
      }
    }

    const std::string & PathManager::GetInputFile() const
    {
      return inputFile;
    }
    const std::string & PathManager::GetImageDirectory() const
    {
      return imageDirectory;
    }
    const std::string & PathManager::GetColloidPath() const
    {
      return colloidFile;
    }
    const std::string & PathManager::GetReportPath() const
    {
      return reportName;
    }

    void PathManager::EmptyOutputDirectories() const
    {
      hemelb::util::DeleteDirContents(imageDirectory);
    }

    hemelb::io::writers::Writer * PathManager::XdrImageWriter(const long int time) const
    {
      char filename[255];
      snprintf(filename, 255, "%08li.dat", time);
#ifdef HEMELB_IMAGES_TO_NULL
      return (new hemelb::io::writers::null::NullWriter());
#else
      return (new hemelb::io::writers::xdr::XdrFileWriter(imageDirectory + std::string(filename)));
#endif
    }

    const std::string& PathManager::GetDataExtractionPath() const
    {
      return dataPath;
    }

    void PathManager::SaveConfiguration(configuration::SimConfig * const simConfig) const
    {
      if (doIo)
      {
        //simConfig->Save(outputDir + "/" + configLeafName);
      }
    }

    void PathManager::GuessOutputDir()
    {
      unsigned long lLastForwardSlash = inputFile.rfind('/');
      if (lLastForwardSlash == std::string::npos)
      {
        // input file supplied is in current folder
        configLeafName = inputFile;
        if (outputDir.length() == 0)
        {
          // no output dir given, defaulting to local.
          outputDir = "./results";
        }
      }
      else
      {
        // input file supplied is a path to the input file
        configLeafName = inputFile.substr(lLastForwardSlash);
        if (outputDir.length() == 0)
        {
          // no output dir given, defaulting to location of input file.
          // note substr is end-exclusive and start-inclusive
          outputDir = inputFile.substr(0, lLastForwardSlash + 1) + "results";
        }
      }
    }
  }
}

