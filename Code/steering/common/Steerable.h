// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_STEERING_COMMON_STEERABLE_H
#define HEMELB_STEERING_COMMON_STEERABLE_H

namespace hemelb
{
  namespace steering
  {
    namespace _private
    {
      // Base class for Steerable member vars
      class SteerableBase
      {
      };
    }

    /**
     * Template class for wrappers around member vars that allows them to
     * be altered by the Steerer. Template parameter should be a suitable
     * subclass of Tag<>; see example in Tags.cc.
     */

    template<class T>
    class Steerable : public _private::SteerableBase
    {
      public:
        typedef T TagClass;
        typedef typename TagClass::WrappedType WrappedType;

      protected:
        // The value we need to steer.
        WrappedType mValue;
        // Constant pointer to a mutable TagClass instance.
        static TagClass* const tag;

      public:
        Steerable() :
            mValue(TagClass::InitialValue)
        {
          Steerable::tag->AddInstance(this);
        }

        ~Steerable()
        {
          tag->RemoveInstance(this);
        }

        WrappedType const& Get(void)
        {
          return mValue;
        }

      protected:
        void Set(WrappedType const& value)
        {
          mValue = value;
        }

        // It's not the Tag that needs to be our friend, it's the Tag's
        // special base class. Happily this is allowed by the standard but
        // friend class TagClass is not.
        friend class TagClass::Super;
    };

    template<class T> T* const Steerable<T>::tag = T::Init();

  }
}
#endif // HEMELB_STEERING_COMMON_STEERABLE_H
