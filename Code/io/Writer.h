#ifndef HEMELB_IO_WRITER_H
#define HEMELB_IO_WRITER_H

#include "lb/LbmParameters.h"
#include "vis/ColPixel.h"
#include "vis/ColourPalette.h"

namespace hemelb
{
  namespace io
  {
    class Writer
    {
      public:

        enum Separator
        {
          eol
        };

        // Special version for eol, using function overloading
        Writer& operator<<(enum Separator const & value);

        // Overload << to write basic types and any necessary separators
        template<typename T>
        Writer& operator<<(T const & value)
        {
          _write(value);
          writeFieldSeparator();
          return *this;
        }

        void writePixel(vis::ColPixel *col_pixel_p,
                        vis::ColourPaletteFunction*colourPalette,
                        vis::DomainStats* iDomainStats,
                        int mode,
                        const lb::StressTypes iLbmStressType);

        // Function to get the current position of writing in the stream.
        virtual unsigned int getCurrentStreamPosition() const = 0;

      protected:
        Writer();
        virtual ~Writer() =0;

        // Functions for formatting control
        virtual void writeFieldSeparator() = 0;
        virtual void writeRecordSeparator() = 0;

        // Methods to simply write (no separators) which are virtual and
        // hence must be overriden.
        virtual void _write(int const & value) = 0;
        virtual void _write(double const & value) = 0;
        virtual void _write(short const & value) = 0;
        virtual void _write(float const & value) = 0;
        virtual void _write(unsigned int const & value) = 0;

    };

  /*template <>
   inline Writer& Writer::operator<< <enum Writer::Separator> (enum Writer::Separator const & value) {
   writeRecordSeparator();
   return *this;
   }*/

  }
}

#endif // HEMELB_IO_WRITER_H
