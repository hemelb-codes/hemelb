// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "reporting/Reporter.h"
#include "ctemplate/template.h"

namespace fs = std::filesystem;

namespace hemelb
{
  namespace reporting
  {
    Reporter::Reporter(const fs::path& rd, const std::string &inputFile) :
        report_dir(rd), dictionary("Reporting dictionary")
    {
      dictionary.SetValue("CONFIG", inputFile);
    }

    void Reporter::AddReportable(Reportable* reportable)
    {
      reportableObjects.push_back(reportable);
    }

    void Reporter::Write(const fs::path& ctemplate, const std::string &as)
    {
      std::string output;
      ctemplate::ExpandTemplate(ctemplate.string(), ctemplate::STRIP_BLANK_LINES, dictionary.GetRaw(), &output);
      auto to = report_dir / as;
      std::fstream file(to, std::ios_base::out);
      file << output << std::flush;
      file.close();
    }

    void Reporter::FillDictionary()
    {
      for (std::vector<Reportable*>::iterator reporters = reportableObjects.begin();
          reporters != reportableObjects.end(); reporters++)
      {
        (*reporters)->Report(dictionary);
      }
    }
  }
}
