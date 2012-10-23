// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_STEERING_COMMON_STEERER_H
#define HEMELB_STEERING_COMMON_STEERER_H

#include <iostream>
#include <map>
#include <string>

#include "steering/common/Singleton.h"
#include "steering/common/Tag.h"

/*
 * Global object to control steerable parameters.
 * 
 * To work, this must have the following methods:
 *
 * SteeringType* Init(void);
 *
 * Returns the steerer object, creating it
 * if need be. Currently do this by inheriting from the Singleton
 * class template.
 *
 * template <typename TagType>
 * void SetTag(typename TagType::WrappedType value)
 *
 * Get the tag type from the list stored and sets its value to the one
 * given.
 *
 * template <typename TagType>
 * void RegisterTag(TagType* tag) 
 *
 * Called by the Tag<> instance to register itself. Must be accessible
 * from there, i.e. make that a friend. However the standard doesn't
 * allow partial specializations to be declared friend so we have to
 * let any Tag<> be our friend, with:
 *
 * template<class T, class DerivedTag, class SteererClass> friend class Tag;
 *
 */
namespace hemelb
{
  namespace steering
  {
    namespace _private
    {
      class TagBase;

      template<class T, class DerivedTag, class SteererClass>
      class Tag;
    }
    class Steerer : public _private::Singleton<Steerer>
    {
      public:
        // Singleton exposes a method static Steerer* Init(void);

        typedef std::map<std::string, _private::TagBase*> ContainerType;
        typedef ContainerType::value_type PairType;

        template<typename TagType>
        void SetTag(typename TagType::WrappedType value)
        {
          TagType* tag = static_cast<TagType*>(mTags[TagType::TagString]);
          tag->SetInstanceValues(value);
        }

      protected:
        template<typename TagType>
        void RegisterTag(TagType* tag)
        {
          mTags.insert(PairType(TagType::TagString, static_cast<_private::TagBase*>(tag)));
        }

        ContainerType mTags;
        Steerer();

        // Singleton must be a friend.
        friend class _private::Singleton<Steerer>;
        // Standard apparently doesn't allow partial specializations to be
        // declared friend so we have to let any Tag<> be our friend.
        template<class T, class DerivedTag, class SteererClass> friend class _private::Tag;
    };

  }
}
#endif // HEMELB_STEERING_COMMON_STEERER_H
