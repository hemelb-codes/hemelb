// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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

    FloatingType LbDataSourceIterator::GetPressure() const
    {
      return converter.ConvertPressureToPhysicalUnits(propertyCache.densityCache.Get(position) * Cs2);
    }

    util::Vector3D<FloatingType> LbDataSourceIterator::GetVelocity() const
    {
      return converter.ConvertVelocityToPhysicalUnits(propertyCache.velocityCache.Get(position));
    }

    FloatingType LbDataSourceIterator::GetShearStress() const
    {
      return converter.ConvertStressToPhysicalUnits(propertyCache.wallShearStressMagnitudeCache.Get(position));
    }

    FloatingType LbDataSourceIterator::GetVonMisesStress() const
    {
      return converter.ConvertStressToPhysicalUnits(propertyCache.vonMisesStressCache.Get(position));
    }

    FloatingType LbDataSourceIterator::GetShearRate() const
    {
      return converter.ConvertShearRateToPhysicalUnits(propertyCache.shearRateCache.Get(position));
    }

    util::Matrix3D LbDataSourceIterator::GetStressTensor() const
    {
      return converter.ConvertFullStressTensorToPhysicalUnits(propertyCache.stressTensorCache.Get(position));
    }

    util::Vector3D<PhysicalStress> LbDataSourceIterator::GetTraction() const
    {
      return converter.ConvertTractionToPhysicalUnits(propertyCache.tractionCache.Get(position),
                                                      data.GetSite(position).GetWallNormal());
    }

    util::Vector3D<PhysicalStress> LbDataSourceIterator::GetTangentialProjectionTraction() const
    {
      return converter.ConvertStressToPhysicalUnits(propertyCache.tangentialProjectionTractionCache.Get(position));
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
      return data.GetProcIdFromGlobalCoords(location) == net::NetworkTopology::Instance()->GetLocalRank();
    }

    distribn_t LbDataSourceIterator::GetVoxelSize() const
    {
      return converter.GetVoxelSize();
    }

    const util::Vector3D<distribn_t>& LbDataSourceIterator::GetOrigin() const
    {
      return converter.GetLatticeOrigin();
    }

    bool LbDataSourceIterator::IsWallSite(const util::Vector3D<site_t>& location) const
    {
      site_t localSiteId = data.GetContiguousSiteId(location);

      return data.GetSite(localSiteId).IsWall();
    }
  }
}
