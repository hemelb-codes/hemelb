#ifndef HEMELB_VIS_SCREENPIXELS_H
#define HEMELB_VIS_SCREENPIXELS_H

#include "io/Writer.h"
#include "ColPixel.h"

namespace hemelb
{
  namespace vis
  {

    struct ScreenPixels
    {
      public:
        ScreenPixels();
        ~ScreenPixels();

        void Reset();

        void FoldIn(const ScreenPixels* inScreen, const VisSettings* visSettings);
        void AddPixel(const ColPixel* newPixel, const VisSettings* visSettings);
        void AddPixels(const ColPixel* newPixel,
                       unsigned int pixelCount,
                       const VisSettings* visSettings);

        void RenderLine(const float endPoint1[3],
                        const float endPoint2[3],
                        const VisSettings* visSettings);

        static const unsigned int COLOURED_PIXELS_MAX = 2048 * 2048;

        void SetSize(int x, int y);

        int GetPixelsX() const;
        int GetPixelsY() const;

        void WritePixels(io::Writer* writer, const DomainStats* domainStats, const VisSettings* visSettings) const;
        void WriteImage(io::Writer* writer,
                        const DomainStats* domainStats,
                        const VisSettings* visSettings) const;

        ColPixel pixels[COLOURED_PIXELS_MAX];
        int pixelId[COLOURED_PIXELS_MAX];
        unsigned int pixelCount;

      private:
        template<bool xLimited> void RenderLineHelper(int x,
                                                      int y,
                                                      int incE,
                                                      int incNE,
                                                      int limit,
                                                      const VisSettings* visSettings);

        unsigned int PixelsX, PixelsY;
    };

  }
}

#endif /* HEMELB_VIS_SCREENPIXELS_H */
