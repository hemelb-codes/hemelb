#include "extraction/LineGeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    LineGeometrySelector::LineGeometrySelector(const util::Vector3D<float>& endpoint1,
                                               const util::Vector3D<float>& endpoint2) :
      endpoint1(endpoint1), lineVector(endpoint2 - endpoint1), lineLength(lineVector.GetMagnitude())
    {

    }

    util::Vector3D<float> LineGeometrySelector::GetEndpoint1() const
    {
      return endpoint1;
    }

    util::Vector3D<float> LineGeometrySelector::GetEndpoint2() const
    {
      return lineVector + endpoint1;
    }

    bool LineGeometrySelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                const util::Vector3D<float>& location)
    {
      const float lengthAlongLine = (lineVector.Dot(location - endpoint1)) / lineLength;

      if (lengthAlongLine < 0. || lengthAlongLine > lineLength)
      {
        return false;
      }

      // Use magnitude squared as it saves a sqrt operation.
      const float perpendicularDistanceSquared =
          ( (endpoint1 + lineVector * lengthAlongLine) - location).GetMagnitudeSquared();

      return perpendicularDistanceSquared <= (0.5 * 0.5 * data.GetVoxelSize() * data.GetVoxelSize());
    }
  }
}
