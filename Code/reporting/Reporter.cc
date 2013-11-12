// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
