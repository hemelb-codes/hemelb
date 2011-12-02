#ifndef HEMELB_CONFIGURATION_COMMANDLINE_H
#define HEMELB_CONFIGURATION_COMMANDLINE_H

#include <string>
#include <stdlib.h>
#include "topology/NetworkTopology.h"
#include "log/Logger.h"
namespace hemelb{
  namespace configuration{
    class CommandLine {
      public:
        CommandLine(int aargc, const char *const *const aargv);
        void PrintUsage();
        unsigned int NumberOfSnapshotsPerCycle() const {
          return(snapshotsPerCycle);
        }
        unsigned int NumberOfImagesPerCycle() const {
          return(imagesPerCycle);
        }
        std::string const & GetOutputDir() const {
           return(outputDir);
        }
        std::string const & GetInputFile() const {
           return(inputFile);
        }
        int GetSteeringSessionId() const {
           return(steeringSessionId);
        }
        int ArgumentCount() const {
          return(argc);
        }
        const char *const * Arguments() {
          return(argv);
        }
        bool HasProblems() {
          return(!ok);
        }
      private:
        std::string inputFile;
        std::string outputDir;
        unsigned int snapshotsPerCycle;
        unsigned int imagesPerCycle;
        int steeringSessionId;
        int argc;
        const char *const *const argv;
        bool ok;
    };
  }
}

#endif //HEMELB_CONFIGURATION_COMMANDLINE_H
