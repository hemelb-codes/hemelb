// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_MULTISCALE_SHAREDVALUE_H
#define HEMELB_MULTISCALE_SHAREDVALUE_H
#include "multiscale/Intercommunicand.h"

namespace hemelb
{
  namespace multiscale
  {
    /***
     * Shared value common base, allowing shared values with different runtime types to be contained together.
     * These will be static cast to specialised shared values, depending on the type of value contained, using runtime type information
     * stored in an IntercommunicandType
     */
    class BaseSharedValue
    {
      public:
        /***
         * Syntactic sugar making it prettier to register shared values with their owner in an initialiser list.
         * @param owner
         */
        BaseSharedValue(Intercommunicand *owner){
          owner->RegisterSharedValue(this);
        }
        
    };
    /***
     * A shared value is a field in an object which can be read and set by the multiscale system.
     * Should be used as a field in classes inheriting Intercommunicand.
     * The cast operators are defined so these fields can be treated as interchangeable with their contents.
     * @tparam payload The type of value contained.
     */
    template<class Payload> class SharedValue : public BaseSharedValue {
      public:
      /***
       * Construct a shared value
       * @param owner The owning intercommunicand, used as "this" in an initialiser list in the parent.
       * @param val An initial value for the payload.
       */
        SharedValue(Intercommunicand* owner, Payload val=Payload())
          :BaseSharedValue(owner),contents(val)
          {}

        /***
         * Cast operator to payload value
         */
        operator Payload &(){
          return contents;
        }
        /***
         * Cast operator to const payload value
         */
        operator const Payload &() const {
          return contents;
        }
        /***
         * Assignment operator from payload value
         * @param input the value to assign
         * @return the value assigned.
         */
        SharedValue<Payload> & operator= (const Payload & input)
          {
            contents=input;
            return *this;
          }
      private:
        Payload contents;
    };
    template<class Payload> std::ostream & operator<<(std::ostream & stream, const SharedValue<Payload> & sv){
      return stream << static_cast<Payload>(sv);
    }
  }
}

#endif // HEMELB_MULTISCALE_SHAREDVALUE_H
