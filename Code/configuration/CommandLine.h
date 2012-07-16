// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_CONFIGURATION_COMMANDLINE_H
#define HEMELB_CONFIGURATION_COMMANDLINE_H

#include <string>

#include "topology/NetworkTopology.h"
#include "log/Logger.h"
namespace hemelb
{
  namespace configuration
  {
    /**
     * Abstraction of HemeLB command line parameters.
     * Arguments should be:
     * - -in input xml configuration file (default input.xml)
     * - -out output folder (empty default, but the hemelb::io::PathManager will guess a value from the input file if not given.)
     * - -i number of images (default 10)
     * - -s number of snapshots (default 10)
     * - -ss steering session i.d. (default 1)
     */
    class CommandLine
    {
      public:
        /**
         * Constructor, constructed from the command line arguments.
         * Parses the command line arguments, and also initialises MPI.
         * @param aargc count of arguments supplied to program, including program name
         * @param aargv values of arguments supplied to program, of which first is program name
         */
        CommandLine(int aargc, const char * const * const aargv);
        /**
         * Report to standard output an error message describing the usage
         */
        void PrintUsage();
        /**
         * @return Number of snapshots that should be produced
         */
        unsigned int NumberOfSnapshots() const
        {
          return (snapshots);
        }
        /**
         * @return Number of images that should be produced
         */
        unsigned int NumberOfImages() const
        {
          return (images);
        }
        /**
         * The output directory that should be used for result files.
         * Empty default, but the hemelb::io::PathManager will guess a value from the input file if not given.)
         * @return Reference to member, the relative or full path to which output files should be written.
         */
        std::string const & GetOutputDir() const
        {
          return (outputDir);
        }
        /**
         * @return Reference to member, the relative or full path from which an input xml config file should be loaded.
         */
        std::string const & GetInputFile() const
        {
          return (inputFile);
        }
        /**
         * @return A unique integer representing the steering session to which to attach.
         */
        int GetSteeringSessionId() const
        {
          return (steeringSessionId);
        }
        /**
         * @return  Total count of command line arguments.
         */
        int ArgumentCount() const
        {
          return (argc);
        }
        /**
         *
         * @return the command line arguments that were given.
         */
        const char * const * Arguments()
        {
          return (argv);
        }
        /**
         * Test this after construction to see if there were problems with creation.
         * Used while we are not using exceptions.
         * @return True if there were problems parsing the command line, false otherwise.
         */
        //TODO replace with an exception
        bool HasProblems()
        {
          return (!ok);
        }
      private:
        std::string inputFile; //! local or full path to input file
        std::string outputDir; //! local or full path to input file
        unsigned int snapshots; //! snapshots to produce
        unsigned int images; //! images to produce
        int steeringSessionId; //! unique identifier for steering session
        int argc; //! count of command line arguments, including program name
        const char * const * const argv; //! command line arguments
        bool ok; //! track if the construction went OK.
    };
  }
}

#endif //HEMELB_CONFIGURATION_COMMANDLINE_H
