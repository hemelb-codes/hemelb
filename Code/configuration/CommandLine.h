#ifndef HEMELB_COMMANDLINE_H
#define HEMELB_COMMANDLINE_H

#include <string>
#include <stdlib.h>
#include "topology/NetworkTopology.h"
#include "log/Logger.h"
namespace hemelb{
  namespace configuration{
    class CommandLine {
      public:
        CommandLine(int aargc, char **aargv);
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
        char ** Arguments() { // Should try to make this const
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
        char ** argv;
        bool ok;
    };
  }
}

#endif
