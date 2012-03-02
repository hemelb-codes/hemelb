#include "extraction/LbDataSourceIterator.h"

namespace hemelb
{
  namespace extraction
  {
    LbDataSourceIterator::LbDataSourceIterator(const lb::MacroscopicPropertyCache& propertyCache,
                                               const geometry::LatticeData& data,
                                               const util::UnitConverter& converter) :
      propertyCache(propertyCache), data(data), converter(converter), position(0)
    {

    }

    bool LbDataSourceIterator::ReadNext(util::Vector3D<site_t>& location,
                                        float& pressure,
                                        util::Vector3D<float>& velocity,
                                        float& stress)
    {
      if (position >= data.GetLocalFluidSiteCount())
      {
        return false;
      }

      ++position;

      location = data.GetSite(position).GetGlobalSiteCoords();
      pressure = converter.ConvertPressureToPhysicalUnits(propertyCache.GetDensity(position) * Cs2);
      velocity = converter.ConvertVelocityToPhysicalUnits(propertyCache.GetVelocity(position));
      stress = converter.ConvertStressToPhysicalUnits(propertyCache.GetStress(position));

      return true;
    }

    void LbDataSourceIterator::Reset()
    {
      position = 0;
    }

    bool LbDataSourceIterator::IsValidLatticeSite(const util::Vector3D<site_t>& location) const
    {
      return data.IsValidLatticeSite(location);
    }

    bool LbDataSourceIterator::IsAvailable(const util::Vector3D<site_t>& location) const
    {
      return data.GetProcIdFromGlobalCoords(location) == topology::NetworkTopology::Instance()->GetLocalRank();
    }

    distribn_t LbDataSourceIterator::GetVoxelSize() const
    {
      return data.GetVoxelSize();
    }
  }
}
