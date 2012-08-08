#ifndef HEMELB_EXTRACTION_OUTPUTFIELD_H
#define HEMELB_EXTRACTION_OUTPUTFIELD_H

namespace hemelb
{
  namespace extraction
  {
    struct OutputField
    {
        enum FieldType
        {
          Pressure,
          Velocity,
          ShearStress,
          VonMisesStress,
          ShearRate,
          StressTensor,
          TractionVector
        };

        std::string name;
        FieldType type;
    };
  }
}

#endif /* HEMELB_EXTRACTION_OUTPUTFIELD_H */
