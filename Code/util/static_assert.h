
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_STATIC_ASSERT_H
#define HEMELB_UTIL_STATIC_ASSERT_H

#ifdef HEMELB_NO_STATIC_ASSERT

#define HEMELB_STATIC_ASSERT(B)

#else // HEMELB_NO_STATIC_ASSERT

#if (__cplusplus > 201103L)

#define HEMELB_STATIC_ASSERT(B) static_assert(B, "Static assertion failed line " __LINE__)

#else

// Fall back to hacky way

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

    // Here we set up compiler specific pragmas to ignore the warning
    // about the typedef we make but do not use.
    
#ifdef __GNUC__
    
    // Note: Clang likes to pretend to be GCC, so this works there too
#define HEMELB_IGNORE_WULT_START					\
    _Pragma("GCC diagnostic push")					\
    _Pragma("GCC diagnostic ignored \"-Wunused-local-typedefs\"")
#define HEMELB_IGNORE_WULT_END			\
    _Pragma("GCC diagnostic pop")

#else
    // Could add further compiler warning supressions here
#define HEMELB_IGNORE_WULT_START
#define HEMELB_IGNORE_WULT_END
#endif
    
#define HEMELB_STATIC_ASSERT( B )					\
    HEMELB_IGNORE_WULT_START						\
    typedef ::hemelb::util::static_assert_test< sizeof(::hemelb::util::STATIC_ASSERTION_FAILURE< HEMELB_UTIL_STATIC_ASSERT_BOOL_CAST( B ) >) > HEMELB_JOIN(static_assert_typedef,__LINE__) \
      HEMELB_IGNORE_WULT_END
  }
}

#endif // HEMELB_CPP11

#endif // HEMELB_NO_STATIC_ASSERT

#endif // HEMELB_UTIL_STATIC_ASSERT_H
