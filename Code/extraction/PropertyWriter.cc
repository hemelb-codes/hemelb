// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "extraction/PropertyWriter.h"

namespace hemelb
{
  namespace extraction
  {
    PropertyWriter::PropertyWriter(IterableDataSource& dataSource,
                                   const std::vector<PropertyOutputFile*>& propertyOutputs)
    {
      for (unsigned outputNumber = 0; outputNumber < propertyOutputs.size(); ++outputNumber)
      {
        localPropertyOutputs.push_back(new LocalPropertyOutput(dataSource, propertyOutputs[outputNumber]));
      }
    }

    PropertyWriter::~PropertyWriter()
    {
      for (unsigned outputNumber = 0; outputNumber < localPropertyOutputs.size(); ++outputNumber)
      {
        delete localPropertyOutputs[outputNumber];
      }
    }

    const std::vector<LocalPropertyOutput*>& PropertyWriter::GetPropertyOutputs() const
    {
      return localPropertyOutputs;
    }

    void PropertyWriter::Write(unsigned long iterationNumber) const
    {
      for (unsigned outputNumber = 0; outputNumber < localPropertyOutputs.size(); ++outputNumber)
      {
        localPropertyOutputs[outputNumber]->Write((uint64_t) iterationNumber);
      }
    }
  }
}
