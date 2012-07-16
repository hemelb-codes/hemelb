// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UTIL_STATIC_ASSERT_H
#define HEMELB_UTIL_STATIC_ASSERT_H

#ifdef HEMELB_NO_STATIC_ASSERT

#define HEMELB_STATIC_ASSERT(B)

#else // HEMELB_NO_STATIC_ASSERT

#define HEMELB_UTIL_STATIC_ASSERT_BOOL_CAST( x ) ((x) == 0 ? false : true)
#define HEMELB_JOIN(X,Y) HEMELB_DO_JOIN(X,Y)
#define HEMELB_DO_JOIN(X,Y) X##Y
namespace hemelb
{
  namespace util
  {

    template<bool x> struct STATIC_ASSERTION_FAILURE;

    template<> struct STATIC_ASSERTION_FAILURE<true>
    {
        enum
        {
          value = 1
        };
    };

    template<int x> struct static_assert_test
    {
    };

#define HEMELB_STATIC_ASSERT( B ) typedef ::hemelb::util::static_assert_test< sizeof(::hemelb::util::STATIC_ASSERTION_FAILURE< HEMELB_UTIL_STATIC_ASSERT_BOOL_CAST( B ) >) > HEMELB_JOIN(static_assert_typedef,__LINE__)
  }
}

#endif // HEMELB_NO_STATIC_ASSERT

#endif // HEMELB_UTIL_STATIC_ASSERT_H
