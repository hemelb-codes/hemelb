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
