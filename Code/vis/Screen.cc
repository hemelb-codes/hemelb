#include <stdlib.h>

#include "topology/NetworkTopology.h"
#include "util/utilityFunctions.h"
#include "vis/Screen.h"

namespace hemelb
{
  namespace vis
  {

    // TODO This is probably going to have to be cleverly redesigned. We need to pass the images around over several iterations without
    // interference between the steering and written-to-disk images.

    Screen::Screen()
    {
    }

    Screen::~Screen()
    {
    }

    /**
     * Add a pixel to the screen.
     *
     * @param newPixel The new pixel to be added
     * @param iStressType The stress type of the visualisation
     * @param mode Controls what aspects of the visualisation to display.
     */
    void Screen::AddPixel(const ColPixel* newPixel, const VisSettings* visSettings)
    {
      pixels.AddPixel(newPixel, visSettings);
    }

    /**
     * Render a line between two points on the screen.
     *
     * @param endPoint1
     * @param endPoint2
     * @param iStressType
     * @param mode
     */
    void Screen::RenderLine(const float endPoint1[3],
                            const float endPoint2[3],
                            const VisSettings* visSettings)
    {
      pixels.RenderLine(endPoint1, endPoint2, visSettings);
    }

    void Screen::Set(float maxX,
                     float maxY,
                     int pixelsX,
                     int pixelsY,
                     float rad,
                     const Viewpoint* viewpoint)
    {
      MaxXValue = maxX;
      MaxYValue = maxX;

      viewpoint->RotateToViewpoint(MaxXValue, 0.0F, 0.0F, UnitVectorProjectionX);
      viewpoint->RotateToViewpoint(0.0F, MaxYValue, 0.0F, UnitVectorProjectionY);

      pixels.SetSize(pixelsX, pixelsY);

      ScaleX = (float) pixels.GetPixelsX() / (2.F * MaxXValue);
      ScaleY = (float) pixels.GetPixelsY() / (2.F * MaxYValue);

      float radVector[3];
      viewpoint->RotateToViewpoint(0.F, 0.F, -rad, radVector);

      for (int ii = 0; ii < 3; ++ii)
      {
        vtx[ii] = (0.5F * radVector[ii]) - UnitVectorProjectionX[ii] - UnitVectorProjectionY[ii];
      }

      UnitVectorProjectionX[0] *= (2.F / (float) pixels.GetPixelsX());
      UnitVectorProjectionX[1] *= (2.F / (float) pixels.GetPixelsX());
      UnitVectorProjectionX[2] *= (2.F / (float) pixels.GetPixelsX());

      UnitVectorProjectionY[0] *= (2.F / (float) pixels.GetPixelsY());
      UnitVectorProjectionY[1] *= (2.F / (float) pixels.GetPixelsY());
      UnitVectorProjectionY[2] *= (2.F / (float) pixels.GetPixelsY());
    }

    void Screen::Resize(unsigned int newPixelsX, unsigned int newPixelsY)
    {
      pixels.SetSize(newPixelsX, newPixelsY);
    }

    void Screen::Reset()
    {
      pixels.pixelCount = 0;
    }

    void Screen::CompositeImage(const VisSettings* visSettings)
    {
      // Status object for MPI comms.
      MPI_Status status;

      /*
       * We do several iterations.
       *
       * On the first, every even proc passes data to the odd proc below, where it is merged.
       * On the second, the difference is two, so proc 3 passes to 1, 7 to 5, 11 to 9 etc.
       * On the third the differenec is four, so proc 5 passes to 1, 13 to 9 etc.
       * .
       * .
       * .
       *
       * This continues until all data is passed back to processor one, which passes it to proc 0.
       */

      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      // Start with a difference in rank of 1, doubling every time.
      for (proc_t deltaRank = 1; deltaRank < netTop->GetProcessorCount(); deltaRank <<= 1)
      {
        // The receiving proc is all the ranks that are 1 modulo (deltaRank * 2)
        for (proc_t receivingProc = 1; receivingProc < (netTop->GetProcessorCount() - deltaRank); receivingProc
            += deltaRank << 1)
        {
          proc_t sendingProc = receivingProc + deltaRank;

          // If we're the sending proc, do the send.
          if (netTop->GetLocalRank() == sendingProc)
          {
            MPI_Send(&pixels.pixelCount,
                     1,
                     MpiDataType(pixels.pixelCount),
                     receivingProc,
                     20,
                     MPI_COMM_WORLD);

            if (pixels.pixelCount > 0)
            {
              MPI_Send(pixels.pixels,
                       pixels.pixelCount,
                       MpiDataType<ColPixel> (),
                       receivingProc,
                       20,
                       MPI_COMM_WORLD);
            }
          }

          // If we're the receiving proc, receive.
          else if (netTop->GetLocalRank() == receivingProc)
          {
            MPI_Recv(&compositingBuffer.pixelCount,
                     1,
                     MpiDataType(compositingBuffer.pixelCount),
                     sendingProc,
                     20,
                     MPI_COMM_WORLD,
                     &status);

            if (compositingBuffer.pixelCount > 0)
            {
              MPI_Recv(compositingBuffer.pixels,
                       compositingBuffer.pixelCount,
                       MpiDataType<ColPixel> (),
                       sendingProc,
                       20,
                       MPI_COMM_WORLD,
                       &status);
            }

            pixels.AddPixels(compositingBuffer.pixels, compositingBuffer.pixelCount, visSettings);
          }
        }
      }

      // Send the final image from proc 1 to 0.
      if (netTop->GetLocalRank() == 1)
      {
        MPI_Send(&pixels.pixelCount, 1, MpiDataType(pixels.pixelCount), 0, 20, MPI_COMM_WORLD);

        if (pixels.pixelCount > 0)
        {
          MPI_Send(pixels.pixels,
                   pixels.pixelCount,
                   MpiDataType<ColPixel> (),
                   0,
                   20,
                   MPI_COMM_WORLD);
        }

      }
      // Receive the final image on proc 0.
      else if (netTop->GetLocalRank() == 0)
      {
        MPI_Recv(&pixels.pixelCount,
                 1,
                 MpiDataType(pixels.pixelCount),
                 1,
                 20,
                 MPI_COMM_WORLD,
                 &status);

        if (pixels.pixelCount > 0)
        {
          MPI_Recv(pixels.pixels,
                   pixels.pixelCount,
                   MpiDataType<ColPixel> (),
                   1,
                   20,
                   MPI_COMM_WORLD,
                   &status);
        }
      }

      pixelCountInBuffer = pixels.pixelCount;

      for (unsigned int m = 0; m < pixelCountInBuffer; m++)
      {
        pixels.pixelId[pixels.pixels[m].GetI() * GetPixelsY() + pixels.pixels[m].GetJ()] = -1;
      }
    }

    bool Screen::MouseIsOverPixel(int mouseX, int mouseY, float* density, float* stress)
    {
      for (unsigned int i = 0; i < pixelCountInBuffer; i++)
      {
        if (pixels.pixels[i].IsRT() && int (pixels.pixels[i].GetI()) == mouseX
            && int (pixels.pixels[i].GetJ()) == mouseY)
        {
          *density = pixels.pixels[i].GetDensity();
          *stress = pixels.pixels[i].GetStress();

          return true;
        }
      }

      return false;
    }

    void Screen::WritePixelCount(io::Writer* writer)
    {
      writer->operator <<(GetPixelsX());
      writer->operator <<(GetPixelsY());
      writer->operator <<(pixelCountInBuffer);
    }

    void Screen::WritePixels(const DomainStats* domainStats,
                             const VisSettings* visSettings,
                             io::Writer* writer)
    {
      int index;
      int pix_data[3];
      unsigned char rgb_data[12];
      int bits_per_char = sizeof(char) * 8;

      for (unsigned int i = 0; i < pixelCountInBuffer; i++)
      {
        ColPixel& col_pixel = pixels.pixels[i];
        // Use a ray-tracer function to get the necessary pixel data.
        col_pixel.rawWritePixel(&index, rgb_data, domainStats, visSettings);

        *writer << index;

        pix_data[0] = (rgb_data[0] << (3 * bits_per_char)) + (rgb_data[1] << (2 * bits_per_char))
            + (rgb_data[2] << bits_per_char) + rgb_data[3];

        pix_data[1] = (rgb_data[4] << (3 * bits_per_char)) + (rgb_data[5] << (2 * bits_per_char))
            + (rgb_data[6] << bits_per_char) + rgb_data[7];

        pix_data[2] = (rgb_data[8] << (3 * bits_per_char)) + (rgb_data[9] << (2 * bits_per_char))
            + (rgb_data[10] << bits_per_char) + rgb_data[11];

        for (int i = 0; i < 3; i++)
        {
          *writer << pix_data[i];
        }
        *writer << io::Writer::eol;
      }
    }

    unsigned int Screen::GetPixelCount() const
    {
      return pixelCountInBuffer;
    }

    const float* Screen::GetVtx() const
    {
      return vtx;
    }
    const float* Screen::GetUnitVectorProjectionX() const
    {
      return UnitVectorProjectionX;
    }
    const float* Screen::GetUnitVectorProjectionY() const
    {
      return UnitVectorProjectionY;
    }
    int Screen::GetPixelsX() const
    {
      return pixels.GetPixelsX();
    }
    int Screen::GetPixelsY() const
    {
      return pixels.GetPixelsY();
    }
  }
}
