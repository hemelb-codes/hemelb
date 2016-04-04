
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
        virtual ~Reportable() {}
        virtual void Report(ctemplate::TemplateDictionary& dictionary) = 0;
    };
  }
}


#endif /* HEMELB_REPORTING_REPORTABLE_H */
