// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXCEPTION_H
#define HEMELB_EXCEPTION_H

#include <sstream>
#include <stdexcept>
#include <string>

namespace hemelb
{
  /**
   * Exception class for HemeLB
   *
   * Implements an ostream interface so you can perform formatted output to
   * it easily like:
   *
   * throw Exception() << "Simulation unstable at t = " << currentTime;
   *
   */
  class Exception : public std::exception
  {
    public:
      /**
       * ctor
       */
      Exception() = default;

      /**
       * Copy ctor required -- note that it appends the old content of the
       * stream to the message of the new one and gives the new one an empty
       * stream.
       *
       * @param that
       */
      Exception(const Exception& that) noexcept : mWhat(that.mStream.str())
      {
      }

      /**
       * Return combined error string.
       * @return
       */
      inline const char *what() const noexcept override
      {
          // TODO: clang's libc++ doesn't support std::ostringstream::view (P0804), so have to alloc to get the data :(
          std::ostringstream tmp;
          std::swap(tmp, mStream);
          mWhat += tmp.str();
          return mWhat.c_str();
      }

      /**
       * Delegate output operations to the Stream.
       * @param t
       * @return
       */
      template<typename T>
      Exception& operator<<(const T& t)
      {
        mStream << t;
        return *this;
      }

    private:
      // Note that these are mutable
      mutable std::ostringstream mStream;
      mutable std::string mWhat;
  };
}
#endif // HEMELB_EXCEPTION_H
