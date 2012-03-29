#include "extraction/LbDataSourceIterator.h"

namespace hemelb
{
  namespace extraction
  {
    LbDataSourceIterator::LbDataSourceIterator(const lb::MacroscopicPropertyCache& propertyCache,
                                               const geometry::LatticeData& data,
                                               const util::UnitConverter& converter) :
        propertyCache(propertyCache), data(data), converter(converter), position(-1)
    {

    }

    bool LbDataSourceIterator::ReadNext()
    {
      ++position;

      if (position >= data.GetLocalFluidSiteCount())
      {
        return false;
      }

      return true;
    }

    util::Vector3D<site_t> LbDataSourceIterator::GetPosition() const
    {
      return data.GetSite(position).GetGlobalSiteCoords();
    }

    float LbDataSourceIterator::GetPressure() const
    {
      return converter.ConvertPressureToPhysicalUnits(propertyCache.densityCache.Get(position) * Cs2);
    }

    util::Vector3D<float> LbDataSourceIterator::GetVelocity() const
    {
      return converter.ConvertVelocityToPhysicalUnits(propertyCache.velocityCache.Get(position));
    }

    float LbDataSourceIterator::GetShearStress() const
    {
      return converter.ConvertStressToPhysicalUnits(propertyCache.shearStressCache.Get(position));
    }

    float LbDataSourceIterator::GetVonMisesStress() const
    {
      return converter.ConvertStressToPhysicalUnits(propertyCache.vonMisesStressCache.Get(position));
    }

    float LbDataSourceIterator::GetShearRate() const
    {
      return converter.ConvertShearRateToPhysicalUnits(propertyCache.shearRateCache.Get(position));
    }

    void LbDataSourceIterator::Reset()
    {
      position = -1;
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
