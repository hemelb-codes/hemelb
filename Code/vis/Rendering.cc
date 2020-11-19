// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "vis/Rendering.h"

namespace hemelb
{
  namespace vis
  {
    Rendering::Rendering(PixelSet<BasicPixel>* glyph, PixelSet<raytracer::RayDataNormal>* ray,
                         PixelSet<streaklinedrawer::StreakPixel>* streak) :
        glyphResult(glyph), rayResult(ray), streakResult(streak)
    {

    }

    void Rendering::ReleaseAll()
    {
      if (glyphResult != nullptr)
      {
        glyphResult->Release();
      }

      if (rayResult != nullptr)
      {
        rayResult->Release();
      }

      if (streakResult != nullptr)
      {
        streakResult->Release();
      }
    }

    void Rendering::ReceivePixelCounts(net::Net* inNet, proc_t source)
    {
      if (glyphResult != nullptr)
      {
        glyphResult->ReceiveQuantity(inNet, source);
      }
      if (rayResult != nullptr)
      {
        rayResult->ReceiveQuantity(inNet, source);
      }
      if (streakResult != nullptr)
      {
        streakResult->ReceiveQuantity(inNet, source);
      }
    }

    void Rendering::ReceivePixelData(net::Net* inNet, proc_t source)
    {
      if (glyphResult != nullptr)
      {
        glyphResult->ReceivePixels(inNet, source);
      }
      if (rayResult != nullptr)
      {
        rayResult->ReceivePixels(inNet, source);
      }
      if (streakResult != nullptr)
      {
        streakResult->ReceivePixels(inNet, source);
      }
    }

    void Rendering::SendPixelCounts(net::Net* inNet, proc_t destination)
    {
      if (glyphResult != nullptr)
      {
        glyphResult->SendQuantity(inNet, destination);
      }
      if (rayResult != nullptr)
      {
        rayResult->SendQuantity(inNet, destination);
      }
      if (streakResult != nullptr)
      {
        streakResult->SendQuantity(inNet, destination);
      }
    }

    void Rendering::SendPixelData(net::Net* inNet, proc_t destination)
    {
      if (glyphResult != nullptr)
      {
        glyphResult->SendPixels(inNet, destination);
      }
      if (rayResult != nullptr)
      {
        rayResult->SendPixels(inNet, destination);
      }
      if (streakResult != nullptr)
      {
        streakResult->SendPixels(inNet, destination);
      }
    }

    void Rendering::Combine(const Rendering& other)
    {
      if (glyphResult != nullptr)
      {
        glyphResult->Combine(*other.glyphResult);
      }
      if (rayResult != nullptr)
      {
        rayResult->Combine(*other.rayResult);
      }
      if (streakResult != nullptr)
      {
        streakResult->Combine(*other.streakResult);
      }
    }

    void Rendering::PopulateResultSet(PixelSet<ResultPixel>* resultSet)
    {
      if (glyphResult != nullptr)
      {
        AddPixelsToResultSet(resultSet, glyphResult->GetPixels());
      }

      if (rayResult != nullptr)
      {
        AddPixelsToResultSet(resultSet, rayResult->GetPixels());
      }

      if (streakResult != nullptr)
      {
        AddPixelsToResultSet(resultSet, streakResult->GetPixels());
      }
    }
  }
}
