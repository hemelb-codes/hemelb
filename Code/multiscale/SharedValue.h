#ifndef HEMELB_MULTISCALE_SHAREDVALUE_H
#define HEMELB_MULTISCALE_SHAREDVALUE_H
#include "multiscale/Intercommunicand.h"

namespace hemelb
{
  namespace multiscale
  {
    class BaseSharedValue
    {
      public:
        BaseSharedValue(Intercommunicand *owner){
          owner->RegisterSharedValue(this);
        }
        
    };
    template<class payload> class SharedValue : public BaseSharedValue {
      public:
        SharedValue(Intercommunicand* owner, payload val=payload())
          :BaseSharedValue(owner),contents(val)
          {}
        payload contents;
        // cast operator
        operator payload &(){
          return contents;
        }
        // assignment operator
        SharedValue<payload> & operator= (const payload & input)
          {
            contents=input;
            return *this;
          }
        
    };
    template<class payload> std::ostream & operator<<(std::ostream & stream, const SharedValue<payload> & sv){
      return stream << sv.contents;
    }
  }
}

#endif // HEMELB_MULTISCALE_SHAREDVALUE_H
