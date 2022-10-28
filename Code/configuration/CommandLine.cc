// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "configuration/CommandLine.h"

namespace hemelb::configuration
{
    CommandLine::CommandLine(int aargc, const char * const aargv[])
            : CommandLine(std::vector<std::string>(aargv, aargv + aargc))
    {
    }

    CommandLine::CommandLine(std::initializer_list<char const*> init)
            : CommandLine(std::vector<std::string>(init.begin(), init.end()))
    {
    }

    CommandLine::CommandLine(std::vector<std::string> const & argv):
            argv(argv)
    {
      // There should be an odd number of arguments since the parameters occur in pairs.
      if ( (argv.size() % 2) == 0)
      {
        throw (OptionError()
            << "There should be an odd number of arguments since the parameters occur in pairs.");
      }

      // All arguments are parsed in pairs, one is a "-<paramName>" type, and one
      // is the <parametervalue>.
      for (size_t ii = 1; ii < argv.size(); ii += 2)
      {
        auto& paramName = argv[ii];
        auto& paramValue = argv[ii + 1];
        if (paramName == "-in")
        {
          inputFile = paramValue;
        }
        else if (paramName == "-out")
        {
          outputDir = paramValue;
        }
        else if (paramName == "-debug")
        {
            if (paramValue == "0") {
                debugMode = false;
            } else if (paramValue == "1") {
                debugMode = true;
            } else {
                throw (OptionError() << "Invalid flag value for -debug");
            }
        }
        else
        {
          throw OptionError() << "Unknown option: " << paramName;
        }
      }

      if (inputFile.empty())
          throw OptionError() << "input file not supplied";

      if (outputDir.empty())
          outputDir = inputFile.parent_path() / "results";

    }

    std::string CommandLine::GetUsage()
    {
        return "Correct usage: hemelb [-<Parameter Name> <Parameter Value>]* \n"
               "Parameter name and significance:\n"
               "\t-in\tPath to the configuration xml file (required)\n"
               "\t-out\tPath to the output folder (default is 'results' in same directory as the input file)\n"
               "\t-debug\tFlag (0 or 1) to enable the hemelb debugger (default: 0)\n";
    }

}
