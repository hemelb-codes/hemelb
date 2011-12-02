#ifndef HEMELB_IO_PATHMANAGER_H
#define HEMELB_IO_PATHMANAGER_H

#include <string>
#include "configuration/CommandLine.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "configuration/SimConfig.h"
#include "io/writers/xdr/XdrFileWriter.h"
namespace hemelb
{
  namespace io
  {
    class PathManager
    {
      public:
        PathManager(const configuration::CommandLine & commandLine,
                    const bool &io,
                    const int &processorCount);
        bool HasProblems() const
        {
          return (!ok);
        }
        const std::string & GetInputFile() const;
        const std::string & GetSnapshotDirectory() const;
        const std::string & GetImageDirectory() const;
        const std::string & GetReportPath() const;
        void SaveConfiguration(configuration::SimConfig * const simConfig) const;
        void EmptyOutputDirectories() const;
        hemelb::io::writers::xdr::XdrFileWriter * XdrImageWriter(const long int time) const;
        const std::string SnapshotPath(unsigned long time) const;
      private:
        void GuessOutputDir();
        void InitialiseReport(const int & processorCount);
        std::string outputDir;
        std::string inputFile;
        std::string snapshotDirectory;
        std::string imageDirectory;
        std::string configLeafName;
        std::string reportName;
        const configuration::CommandLine &options;
        bool ok;
        bool doIo;
    };
  }
}

#endif //HEMELB_IO_PATHMANAGER_H
