#ifndef HEMELB_VIS_VIEWPOINT_H
#define HEMELB_VIS_VIEWPOINT_H

namespace hemelb
{
  namespace vis
  {
    class Viewpoint
    {
      public:
        void RotateToViewpoint(float iXIn, float iYIn, float iZIn, float rotatedVector[3]) const;

        void Project(const float p1[], float p2[]) const;

        void SetViewpointPosition(float longitude,
                                  float latitude,
                                  float localCentre[3],
                                  float rad,
                                  float distance);

        const float* GetViewpointCentre() const;

      private:
        void Rotate(float sinX,
                    float cosX,
                    float sinY,
                    float cosY,
                    float xIn,
                    float yIn,
                    float zIn,
                    float rotatedVector[3]) const;

        float dist;
        float SinYRotation, CosYRotation;
        float SinXRotation, CosXRotation;
        float x[3];
    };
  }
}

#endif /* HEMELB_VIS_VIEWPOINT_H */
