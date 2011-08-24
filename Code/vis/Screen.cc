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
      mPixels= new ScreenPixels();
    }

    Screen::~Screen()
    {
      delete mPixels;
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
      mPixels->AddPixel(newPixel, visSettings);
    }

    /**
     * Render a line between two points on the screen.
     *
     * @param endPoint1
     * @param endPoint2
     * @param iStressType
     * @param mode
     */
    void Screen::RenderLine(const XYCoordinates<float>& endPoint1,
                            const XYCoordinates<float>& endPoint2,
                            const VisSettings* visSettings)
    {
      mPixels->RenderLine(endPoint1, endPoint2, visSettings);
    }

    void Screen::Set(float maxX,
                     float maxY,
                     int pixelsX,
                     int pixelsY,
                     float rad,
                     const Viewpoint* iViewpoint)
    {
      mMaxXValue = maxX;
      mMaxYValue = maxX;

      mPixels->SetSize(pixelsX, pixelsY);

      mPixelsPerUnitX = (float) mPixels->GetPixelsX() / (2.F * mMaxXValue);
      mPixelsPerUnitY = (float) mPixels->GetPixelsY() / (2.F * mMaxYValue);

      Vector3D<float> lCameraToLocalCentreVector = iViewpoint->
	RotateCameraCoordinatesToWorldCoordinates
	(Vector3D<float>(0.F, 0.F, -iViewpoint->mDistanceFromEyeToScreen));

      Vector3D<float> lMiddleCentreToMiddleRightOfScreen = iViewpoint->
	RotateCameraCoordinatesToWorldCoordinates(Vector3D<float>(mMaxXValue, 0.0F, 0.0F));

      Vector3D<float> lLowerCentreToTopCentreOfScreen = iViewpoint->
	RotateCameraCoordinatesToWorldCoordinates(Vector3D<float>(0.0F, mMaxYValue, 0.0F));

      mCameraToBottomLeftOfScreen = (lCameraToLocalCentreVector - lMiddleCentreToMiddleRightOfScreen) - lLowerCentreToTopCentreOfScreen;

      mPixelUnitVectorProjectionX = lMiddleCentreToMiddleRightOfScreen * (2.F / (float) mPixels->GetPixelsX());
  
      mPixelUnitVectorProjectionY = lLowerCentreToTopCentreOfScreen * (2.F / (float) mPixels->GetPixelsY());
    }

    void Screen::Resize(unsigned int newPixelsX, unsigned int newPixelsY)
    {
      mPixels->SetSize(newPixelsX, newPixelsY);
    }

    void Screen::Reset()
    {
      mPixels->Reset();
    }

    bool Screen::MouseIsOverPixel(int mouseX, int mouseY, float* density, float* stress)
    {
      const ColPixel* screenPix = mPixels->GetPixelArray();

      for (unsigned int i = 0; i < mPixels->GetStoredPixelCount(); i++)
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
      return mPixels->GetStoredPixelCount();
    }

    const Vector3D<float>& Screen::GetCameraToBottomLeftOfScreenVector() const
    {
      return mCameraToBottomLeftOfScreen;
    }
    const Vector3D<float>& Screen::GetPixelUnitVectorProjectionX() const
    {
      return mPixelUnitVectorProjectionX;
    }
    const Vector3D<float>& Screen::GetPixelUnitVectorProjectionY() const
    {
      return mPixelUnitVectorProjectionY;
    }
    int Screen::GetPixelsX() const
    {
      return mPixels->GetPixelsX();
    }
    int Screen::GetPixelsY() const
    {
      return mPixels->GetPixelsY();
    }

    ScreenPixels* Screen::SwapBuffers(ScreenPixels* inPix)
    {
      ScreenPixels* temp = mPixels;
      mPixels = inPix;
      return temp;
    }
    const ScreenPixels* Screen::GetPixels() const
    {
      return mPixels;
    }
  }
}
