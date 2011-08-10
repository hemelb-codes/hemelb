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
      pixels = new ScreenPixels();
    }

    Screen::~Screen()
    {
      delete pixels;
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
      pixels->AddPixel(newPixel, visSettings);
    }

    /**
     * Render a line between two points on the screen.
     *
     * @param endPoint1
     * @param endPoint2
     * @param iStressType
     * @param mode
     */
    void Screen::RenderLine(const Location<float>& endPoint1,
                            const Location<float>& endPoint2,
                            const VisSettings* visSettings)
    {
      pixels->RenderLine(endPoint1, endPoint2, visSettings);
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

      UnitVectorProjectionX = viewpoint->RotateToViewpoint(Location<float>(MaxXValue, 0.0F, 0.0F));
      UnitVectorProjectionY = viewpoint->RotateToViewpoint(Location<float>(0.0F, MaxYValue, 0.0F));

      pixels->SetSize(pixelsX, pixelsY);

      ScaleX = (float) pixels->GetPixelsX() / (2.F * MaxXValue);
      ScaleY = (float) pixels->GetPixelsY() / (2.F * MaxYValue);

      Location<float> radVector = viewpoint->RotateToViewpoint(Location<float>(0.F, 0.F, -rad));

      mVtx = radVector*0.5F - UnitVectorProjectionX - UnitVectorProjectionY;

      UnitVectorProjectionX.x *= (2.F / (float) pixels->GetPixelsX());
      UnitVectorProjectionX.y *= (2.F / (float) pixels->GetPixelsX());
      UnitVectorProjectionX.z *= (2.F / (float) pixels->GetPixelsX());

      UnitVectorProjectionY.x *= (2.F / (float) pixels->GetPixelsY());
      UnitVectorProjectionY.y *= (2.F / (float) pixels->GetPixelsY());
      UnitVectorProjectionY.z *= (2.F / (float) pixels->GetPixelsY());
    }

    void Screen::Resize(unsigned int newPixelsX, unsigned int newPixelsY)
    {
      pixels->SetSize(newPixelsX, newPixelsY);
    }

    void Screen::Reset()
    {
      pixels->Reset();
    }

    bool Screen::MouseIsOverPixel(int mouseX, int mouseY, float* density, float* stress)
    {
      const ColPixel* screenPix = pixels->GetPixelArray();

      for (unsigned int i = 0; i < pixels->GetStoredPixelCount(); i++)
      {
        if (screenPix[i].IsRT() && int (screenPix[i].GetI()) == mouseX && int (screenPix[i].GetJ())
            == mouseY)
        {
          *density = screenPix[i].GetDensity();
          *stress = screenPix[i].GetStress();

          return true;
        }
      }

      return false;
    }

    unsigned int Screen::GetPixelCount() const
    {
      return pixels->GetStoredPixelCount();
    }

    const Location<float>& Screen::GetVtx() const
    {
      return mVtx;
    }
    const Location<float>& Screen::GetUnitVectorProjectionX() const
    {
      return UnitVectorProjectionX;
    }
    const Location<float>& Screen::GetUnitVectorProjectionY() const
    {
      return UnitVectorProjectionY;
    }
    int Screen::GetPixelsX() const
    {
      return pixels->GetPixelsX();
    }
    int Screen::GetPixelsY() const
    {
      return pixels->GetPixelsY();
    }

    ScreenPixels* Screen::SwapBuffers(ScreenPixels* inPix)
    {
      ScreenPixels* temp = pixels;
      pixels = inPix;
      return temp;
    }
    const ScreenPixels* Screen::GetPixels() const
    {
      return pixels;
    }
  }
}
