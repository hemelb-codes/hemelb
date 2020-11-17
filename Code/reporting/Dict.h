// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REPORTING_DICT_H
#define HEMELB_REPORTING_DICT_H

#include <memory>
#include <string>

// Forward declare
namespace ctemplate {
  class TemplateDictionary;
}

namespace hemelb
{
  namespace reporting
  {
    // Minimal wrapper around ctemplate::TemplateDictionary
    // It simply delegates to that instance
    class Dict {
    public:
      Dict(const std::string&);
      Dict AddSectionDictionary(const std::string&);
      
      void SetValue(const std::string& variable, const std::string& value);
      void SetIntValue(const std::string& variable, long value);

      template<typename T>
      void SetFormattedValue(const std::string& variable, const char* format, const T& value);

      const ctemplate::TemplateDictionary* GetRaw() const;
    private:
      // Private constructor for wrapping sub dicts
      Dict(ctemplate::TemplateDictionary*);
      // Deleter type (for expressing ownership of the wrapped instance)
      typedef void (*Deleter)(ctemplate::TemplateDictionary*);
      std::unique_ptr<ctemplate::TemplateDictionary, Deleter> raw;
    };
  }
}
#endif
