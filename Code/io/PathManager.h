// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_PATHMANAGER_H
#define HEMELB_IO_PATHMANAGER_H

#include <string>
#include "util/fileutils.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace configuration
  {
    class CommandLine;
    class SimConfig;
  }

  namespace io
  {
    // Forward declare the Writer
    namespace writers
    {
      class Writer;
    }

    /**
     * Manage the input and output file system locations for HemeLB reports, extracted data, and input xml.
     */
    class PathManager
    {
      public:
        /**
         * During construction, the path manager will guess an appropriate path from the input file path given at the command line,
         * if one was not specified explicitly.
         * @param commandLine The command-line parser object instantiated from the arguments to the C main function.
         * @param io Set to true if this MPI node is the I/O process and files should be written from this node.
         * @param processorCount The total count of processors, used in generating file-names for report files.
         */
        PathManager(const configuration::CommandLine & commandLine, const bool &io,
                    const int &processorCount);

        /**
         * A local or full path to the input xml configuration file.
         * @return A local or full path to the input xml configuration file.
         */
        const std::string & GetInputFile() const;
        /**
         * Gets the path to the file where colloid output should be written
         * @return
         */
        const std::string & GetColloidPath() const;
        /**
         * Path to where a run report file should be created.
         * @return Reference to path to where a run report file should be created.
         */
        const std::string & GetReportPath() const;
        /**
         * Save the current configuration as a configuration xml file to the output folder.
         * @param simConfig The input configuration instance, constructed from the input xml file.
         */
        void SaveConfiguration(configuration::SimConfig * const simConfig) const;
        /**
         * Delete the output data.
         * This deletes config files, and the report.
         * It is used if the simulation needs to reset.
         */
        void EmptyOutputDirectories() const;

        /**
         * Return the path that property extraction output should go to.
         * @return
         */
        const std::string& GetDataExtractionPath() const;

        /**
         * Create a subdirectory inside the RBC output directory and return its path
         * @param subdirectoryName Name of the subdirectory to be created
         * @return Path to the newly created subdirectory
         */
        const std::string GetRBCOutputPathWithSubdir(std::string subdirectoryName) const;

      private:
        void GuessOutputDir(); //! String processing to generate an appropriate outptu folder name.
        std::string outputDir;
        std::string inputFile;
        std::string colloidFile;
        std::string configLeafName;
        std::string reportName;
        std::string dataPath;
        const configuration::CommandLine &options;
        bool doIo; //! Am I the input/output node?
        std::string rbcsPath; //! Path for RBC output
    };
  }
}

#endif //HEMELB_IO_PATHMANAGER_H
