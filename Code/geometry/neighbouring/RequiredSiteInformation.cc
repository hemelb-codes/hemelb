
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/neighbouring/RequiredSiteInformation.h"
#include "net/mpi.h"
namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {
      RequiredSiteInformation::RequiredSiteInformation(bool initial) :
          choices(terms::Length, initial)
      {
      }
      void RequiredSiteInformation::Require(terms::Term term)
      {
        choices[term] = true;
      }
      bool RequiredSiteInformation::RequiresAny()
      {
        for (std::vector<bool>::iterator choice = choices.begin(); choice != choices.end(); choice++)
        {
          if (*choice)
          {
            return true;
          }
        }
        return false;
      }
      bool RequiredSiteInformation::RequiresAnyFieldDependent()
      {
        for (int choice = terms::Distribution; choice < terms::Length; choice++)
        {
          if (choices[choice])
          {
            return true;
          }
        }
        return false;
      }
      bool RequiredSiteInformation::RequiresAnyNonFieldDependent()
      {
        for (int choice = terms::SiteData; choice <= terms::WallNormal; choice++)
        {
          if (choices[choice])
          {
            return true;
          }
        }
        return false;
      }
      bool RequiredSiteInformation::RequiresAnyMacroscopic()
      {
        for (int choice = terms::Velocity; choice < terms::Length; choice++)
        {
          if (choices[choice])
          {
            return true;
          }
        }
        return false;
      }

      void RequiredSiteInformation::Or(const RequiredSiteInformation& other)
      {
        for (int choice = terms::SiteData; choice < terms::Length; choice++)
        {
          choices[choice]=choices[choice] || other.choices[choice];
        }
      }
      void RequiredSiteInformation::And(const RequiredSiteInformation& other)
      {
        for (int choice = terms::SiteData; choice < terms::Length; choice++)
        {
         choices[choice]=choices[choice] && other.choices[choice];
        }
      }
    }
  }
}
