// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "configuration/CommandLine.h"
#include <cstring>
#include <cstdlib>
namespace hemelb
{
  namespace configuration
  {

    CommandLine::CommandLine(std::vector<std::string> const & argv):
        inputFile("input.xml"), outputDir(""), images(10), steeringSessionId(1), debugMode(false),
      argv(argv)
    {
      // There should be an odd number of arguments since the parameters occur in pairs.
      if ( (argv.size() % 2) == 0)
      {
        throw OptionError()
            << "There should be an odd number of arguments since the parameters occur in pairs.";
      }

      // All arguments are parsed in pairs, one is a "-<paramName>" type, and one
      // is the <parametervalue>.
      for (size_t ii = 1; ii < argv.size(); ii += 2)
      {
        auto paramName = argv[ii].c_str();
        auto paramValue = argv[ii + 1].c_str();
        if (std::strcmp(paramName, "-in") == 0)
        {
          inputFile = std::string(paramValue);
        }
        else if (std::strcmp(paramName, "-out") == 0)
        {
          outputDir = std::string(paramValue);
        }
        else if (std::strcmp(paramName, "-i") == 0)
        {
          char *dummy;
          images = (unsigned int) (strtoul(paramValue, &dummy, 10));
        }
        else if (std::strcmp(paramName, "-ss") == 0)
        {
          char *dummy;
          steeringSessionId = (unsigned int) (strtoul(paramValue, &dummy, 10));
        }
        else if (std::strcmp(paramName, "-debug") == 0)
        {
          debugMode = std::strcmp(paramName, "0") == 0 ?
            false :
            true;
        }
        else
        {
          throw OptionError() << "Unknown option: " << paramName;
        }
      }
    }

    std::string CommandLine::GetUsage()
    {
      std::string ans("Correct usage: hemelb [-<Parameter Name> <Parameter Value>]* \n");
      ans.append("Parameter name and significance:\n");
      ans.append("-in \t Path to the configuration xml file (default is config.xml)\n");
      ans.append("-out \t Path to the output folder (default is based on input file, e.g. config_xml_results)\n");
      ans.append("-i \t Number of images to create (default is 10)\n");
      ans.append("-ss \t Steering session identifier (default is 1)\n");
      return ans;
    }
  }
}
