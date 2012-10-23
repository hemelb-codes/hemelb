// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_GEOMETRY_NEIGHBOURING_REQUIREDSITEINFORMATION_H
#define HEMELB_GEOMETRY_NEIGHBOURING_REQUIREDSITEINFORMATION_H
#include <vector>
#include <algorithm>
namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {

      namespace terms
      {
        enum Term
        {
          SiteData = 0,
          WallDistance = 1,
          WallNormal = 2,
          Distribution = 3,
          Velocity = 4,
          Density = 5,
          ShearStress = 6,
          VonMisesStress = 7,
          ShearRate = 8,
          Length
        };
      }
      /***
       * Class to represent, eventually, exactly what information is needed from a remote site.
       * Initially, this will not be used -- we will just transfer all information.
       */
      class RequiredSiteInformation
      {
        public:
          RequiredSiteInformation(bool initial = false);
          void Or(const RequiredSiteInformation& other);
          void And(const RequiredSiteInformation& other);
          void Require(terms::Term term);
          bool RequiresAny();
          bool RequiresAnyFieldDependent();
          bool RequiresAnyNonFieldDependent();
          bool RequiresAnyMacroscopic();
          bool IsRequired(terms::Term term)
          {
            return choices[term];
          }
        private:
          std::vector<bool> choices;
      };
    }
  }
}
#endif
