// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "reporting/Dict.h"
#include <ctemplate/template.h>

namespace hemelb
{
  namespace reporting
  {
    namespace
    {
      void RealDeleter(ctemplate::TemplateDictionary* ptr)
      {
        delete ptr;
      }
      void NullDeleter(ctemplate::TemplateDictionary* ptr)
      {
      }
      
    }
    Dict::Dict(const std::string& s) : raw(new ctemplate::TemplateDictionary(s), RealDeleter)
    {
    }

    Dict::Dict(ctemplate::TemplateDictionary* tmp) : raw(tmp, NullDeleter)
    {
    }

    const ctemplate::TemplateDictionary* Dict::GetRaw() const
    {
      return raw.get();
    }
    
    Dict Dict::AddSectionDictionary(const std::string& s)
    {
      return Dict(raw->AddSectionDictionary(s));
    }

    void Dict::SetValue(const std::string& variable, const std::string& value)
    {
      raw->SetValue(variable, value);
    }

    void Dict::SetIntValue(const std::string& variable, long value)
    {
      raw->SetIntValue(variable, value);
    }
    
    template <typename T>
    void Dict::SetFormattedValue(const std::string& variable, const char* format, const T& value)
    {
      raw->SetFormattedValue(variable, format, value);
    }

    template 
    void Dict::SetFormattedValue<double>(const std::string& variable, const char* format, const double& value);

  }
}
