#ifndef HEMELB_VIS_RAYTRACER_SITEDATA_H
#define HEMELB_VIS_RAYTRACER_SITEDATA_H

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      //Stores the data about an individual voxel 
      struct SiteData_t
      {
        public:
          SiteData_t(float iValue) :
              mDensity(iValue), mVelocity(iValue), mStress(iValue)
          {
          }

          float GetDensity() const
          {
            return mDensity;
          }

          void SetDensity(float iDensity)
          {
            mDensity = iDensity;
          }

          float GetVelocity() const
          {
            return mVelocity;
          }

          void SetVelocity(float iVelocity)
          {
            mVelocity = iVelocity;
          }

          float GetStress() const
          {
            return mStress;
          }

          void SetStress(float iStress)
          {
            mStress = iStress;
          }

        private:
          float mDensity;
          float mVelocity;
          float mStress;
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_SITEDATA_H
