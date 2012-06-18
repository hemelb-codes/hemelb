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
      class RequiredSiteInformation
      {
        public:
          RequiredSiteInformation();
          void Or(const RequiredSiteInformation& other);
          void And(const RequiredSiteInformation& other);
          void Require(terms::Term term);
          bool Any();
          bool AnyFieldDependent();
          bool AnyNonFieldDependent();
          bool AnyMacroscopic();
          bool Required(terms::Term term){return choices[term];}
        private:
          std::vector<bool> choices;
      };
    }
  }
}
#endif
