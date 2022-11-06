// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_CONFIGURATION_COMMANDLINE_H
#define HEMELB_CONFIGURATION_COMMANDLINE_H

#include <vector>
#include <string>
#include <filesystem>

#include "Exception.h"

namespace hemelb::configuration
{
    /**
     * Abstraction of HemeLB command line parameters.
     * Arguments should be:
     * - -in input xml configuration file (required)
     * - -out output folder (default "results")
     */
    class CommandLine
    {
    private:
        std::filesystem::path inputFile; //! local or full path to input file
        std::filesystem::path outputDir; //! local or full path to output directory
        bool debugMode = false; //! Use debugger
        std::vector<std::string> argv; //! command line arguments

    public:
        /**
         * Allow separate catching of these exceptions to enable usage printing.
         */
        class OptionError : public Exception
        {
        };
        /**
         * Report to standard output an error message describing the usage
         */
        static std::string GetUsage();

        /**
         * Constructed from the command line arguments.
         * @param aargc count of arguments supplied to program, including program name
         * @param aargv values of arguments supplied to program, of which first is program name
         */
        CommandLine(int aargc, const char * const aargv[]);

        CommandLine(std::initializer_list<char const*> init);

        explicit CommandLine(std::vector<std::string> const &argv);

        /**
         * The output directory that should be used for result files.
         * Empty default, but the hemelb::io::PathManager will guess a value from the input file if not given.)
         * @return Reference to member, the relative or full path to which output files should be written.
         */
        [[nodiscard]] inline std::filesystem::path const& GetOutputDir() const
        {
            return outputDir;
        }

        /**
         * @return Reference to member, the relative or full path from which an input xml config file should be loaded.
         */
        [[nodiscard]] inline std::filesystem::path const & GetInputFile() const
        {
            return inputFile;
        }

        /**
         * @return Whether the user requested a debug mode.
         */
        [[nodiscard]] inline bool GetDebug() const
        {
            return debugMode;
        }

        /**
         * @return  Total count of command line arguments.
         */
        [[nodiscard]] inline auto ArgumentCount() const
        {
            return argv.size();
        }
        /**
         *
         * @return the command line arguments that were given.
         */
        [[nodiscard]] inline std::vector<std::string> const& Arguments() const
        {
            return argv;
        }
    };
}

#endif //HEMELB_CONFIGURATION_COMMANDLINE_H
