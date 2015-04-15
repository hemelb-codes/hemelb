// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_STEERING_COMMON_TAG_H
#define HEMELB_STEERING_COMMON_TAG_H

#include <string>
#include <set>

#include "steering/common/Singleton.h"
#include "steering/common/TagBase.h"

#include "steering/common/Steerer.h"
#include "steering/common/Steerable.h"

//class Steerer;
namespace hemelb
{
  namespace steering
  {
    namespace _private
    {
      /**
       * Family of Tags for a type of variable T.  Doesn't specify a tag
       * string; that's the job of the the subclass. Hence we're a "static
       * abstract base class" and using the Curiously Recurring Template
       * Pattern/static polymorphism. Note this is without the common
       * specialization when we don't have a subclass.
       */

      template<typename T, typename DerivedTag, typename SteererClass = class Steerer>
      class Tag : public TagBase,
                  public Singleton<Tag<T, DerivedTag, SteererClass> >
      {
        public:
          // Typedefs to avoid lots of parameterized templates
          typedef T WrappedType;
          typedef Tag<WrappedType, DerivedTag, SteererClass> TagType;
          typedef Steerable<DerivedTag> SteerableType;
          typedef std::set<SteerableType*> ContainerType;
          // String for the tag
          static const std::string TagString;
          static const WrappedType InitialValue;

        protected:
          // Contains of pointers to instances of SteerableType
          ContainerType Instances;
          // Our global Steerer.
          static SteererClass* steerer;

          static DerivedTag* Init()
          {
            // Wrap the singleton behaviour, registering the instance with the
            // steerer. Note have to use static polymorphism pattern here.
            DerivedTag* self =
                static_cast<DerivedTag*>(Singleton<Tag<T, DerivedTag, SteererClass> >::Init());

            if (steerer == nullptr)
              steerer = SteererClass::Init();

            steerer->SteererClass::RegisterTag(self);
            return self;
          }

          // Called by constructor of SteerableType
          void AddInstance(SteerableType* const inst)
          {
            Instances.insert(inst);
          }
          // Called by destructor of SteerableType
          void RemoveInstance(SteerableType* const inst)
          {
            Instances.erase(inst);
          }

        public:
          // Called by the Steerer to update our instances with the new value
          void SetInstanceValues(const T value)
          {
            for (typename ContainerType::iterator i = Instances.begin(); i != Instances.end(); ++i)
            {
              (*i)->Set(value);
            }
          }

          friend class Singleton<Tag<T, DerivedTag, SteererClass> > ;
          friend class Steerable<DerivedTag> ;

          // Standard doesn't allow friend class SteererClass.
      };

      // Template static initializers
      template<typename T, typename DerivedTag, typename SteererClass>
      SteererClass* Tag<T, DerivedTag, SteererClass>::steerer = SteererClass::Init();

      template<typename T, typename DerivedTag, typename SteererClass>
      const std::string Tag<T, DerivedTag, SteererClass>::TagString = "";

      template<typename T, typename DerivedTag, typename SteererClass>
      const typename Tag<T, DerivedTag, SteererClass>::WrappedType Tag<T, DerivedTag, SteererClass>::InitialValue =
          0;

    }
  }
}

#define HEMELB_STEERING_DECLARE_TAG(name, type)                                 \
class name :                                                                    \
  public hemelb::steering::_private::Tag<type, name, hemelb::steering::Steerer> \
  {                                                                             \
  public:                                                                       \
    typedef hemelb::steering::_private::Tag<type, name, Steerer> Super;         \
  }

#define HEMELB_STEERING_INITIALISE_TAG(name, value)                             \
  template<> const std::string name::Super::TagString = #name;                  \
  template<> const name::Super::WrappedType name::Super::InitialValue = value

#endif // HEMELB_STEERING_COMMON_TAG_H
