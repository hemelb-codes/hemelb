// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_REPORTING_REPORTABLE_H
#define HEMELB_REPORTING_REPORTABLE_H

#include "ctemplate/template.h"

namespace hemelb
{
  namespace reporting
  {
    /**
     * Defines the interface for classes that can have entries in reports.
     */
    class Reportable
    {
      public:
        virtual void Report(ctemplate::TemplateDictionary& dictionary) = 0;
        virtual ~Reportable()
        {

        }
    };
  }
}


#endif /* HEMELB_REPORTING_REPORTABLE_H */
