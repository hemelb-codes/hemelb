
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_UTILITYSTRUCTS_H
#define HEMELB_UTIL_UTILITYSTRUCTS_H

namespace hemelb
{
  namespace util
  {

    // Allows sorting key-value pairs using the standard library sort
    template<typename keyType, typename valueType>
    struct key_value_pair
    {
      public:
        keyType key;
        valueType value;

        bool operator<(const key_value_pair other_key_value_pair) const
        {
          return key < other_key_value_pair.key;
        }
    };

  }
}

#endif /* HEMELB_UTIL_UTILITYSTRUCTS_H */
