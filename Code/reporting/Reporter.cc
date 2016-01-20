
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "reporting/Reporter.h"

namespace hemelb
{
  namespace reporting
  {
    Reporter::Reporter(const std::string &apath, const std::string &inputFile) :
        path(apath), imageCount(0), dictionary("Reporting dictionary")
    {
      dictionary.SetValue("CONFIG", inputFile);
    }

    void Reporter::Image()
    {
      imageCount++;
    }

    void Reporter::AddReportable(Reportable* reportable)
    {
      reportableObjects.push_back(reportable);
    }

    void Reporter::Write(const std::string &ctemplate, const std::string &as)
    {
      std::string output;
      ctemplate::ExpandTemplate(ctemplate, ctemplate::STRIP_BLANK_LINES, &dictionary, &output);
      std::string to = path + "/" + as;
      std::fstream file(to.c_str(), std::ios_base::out);
      file << output << std::flush;
      file.close();
    }
    ;

    void Reporter::FillDictionary()
    {
      dictionary.SetIntValue("IMAGES", imageCount);

      for (std::vector<Reportable*>::iterator reporters = reportableObjects.begin();
          reporters != reportableObjects.end(); reporters++)
      {
        (*reporters)->Report(dictionary);
      }
    }
  }
}
