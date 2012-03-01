#include "extraction/SnapshotParser.h"

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

    bool SnapshotParser::ReadNext(util::Vector3D<float>& position,
                                  float& pressure,
                                  util::Vector3D<float>& velocity,
                                  float& stress)
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

      position = (coords + minimalSite) * voxelSize + origin;

      reader.readFloat(pressure);
      reader.readFloat(velocity.x);
      reader.readFloat(velocity.y);
      reader.readFloat(velocity.z);

      reader.readFloat(stress);

      return true;
    }
  }
}
