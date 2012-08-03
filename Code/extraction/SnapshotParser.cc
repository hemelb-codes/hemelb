#include "extraction/SnapshotParser.h"
#include <cassert>

namespace hemelb
{
  namespace extraction
  {
    SnapshotParser::SnapshotParser(std::string filename) :
        file(fopen(filename.c_str(), "r")), reader(file), sitesRead(0)
    {
      unsigned hemelbMagicNumber, snapshotMagicNumber, snapshotVersionNumber;
      reader.readUnsignedInt(hemelbMagicNumber);
      reader.readUnsignedInt(snapshotMagicNumber);
      reader.readUnsignedInt(snapshotVersionNumber);

      unsigned headerLength;
      reader.readUnsignedInt(headerLength);

      int stability;
      reader.readInt(stability);

      reader.readDouble(voxelSize);
      reader.readDouble(origin.x);
      reader.readDouble(origin.y);
      reader.readDouble(origin.z);

      reader.readInt(minimalSite.x);
      reader.readInt(minimalSite.y);
      reader.readInt(minimalSite.z);

      util::Vector3D<int> siteMax;
      reader.readInt(siteMax.x);
      reader.readInt(siteMax.y);
      reader.readInt(siteMax.z);

      reader.readInt(fluidSiteCount);
    }

    SnapshotParser::~SnapshotParser()
    {
      if (file != NULL)
      {
        fclose(file);
      }
    }

    bool SnapshotParser::ReadNext()
    {
      if (sitesRead >= fluidSiteCount)
      {
        return false;
      }

      sitesRead++;

      util::Vector3D<int> coords;
      reader.readInt(coords.x);
      reader.readInt(coords.y);
      reader.readInt(coords.z);

      readPosition = (coords + minimalSite);

      reader.readFloat(readPressure);
      reader.readFloat(readVelocity.x);
      reader.readFloat(readVelocity.y);
      reader.readFloat(readVelocity.z);

      reader.readFloat(readShearStress);

      return true;
    }

    util::Vector3D<site_t> SnapshotParser::GetPosition() const
    {
      return readPosition;
    }

    FloatingType SnapshotParser::GetPressure() const
    {
      return readPressure;
    }

    util::Vector3D<FloatingType> SnapshotParser::GetVelocity() const
    {
      return readVelocity;
    }

    FloatingType SnapshotParser::GetShearStress() const
    {
      return readShearStress;
    }

    FloatingType SnapshotParser::GetVonMisesStress() const
    {
      return readVonMisesStress;
    }

    bool SnapshotParser::IsEdgeSite(const util::Vector3D<site_t>& location) const
    {
      // The current implementation of SnapshotParser does not know enough about the domain to be able implement this method
      assert(false);

      return false;
    }

  }
}
