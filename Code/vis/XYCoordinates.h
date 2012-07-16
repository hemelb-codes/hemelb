// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_XYCOORDINATES_H
#define HEMELB_VIS_XYCOORDINATES_H

#include "constants.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace vis
  {
    template<class T>
    class XYCoordinates
    {
      public:
        T x, y;

        XYCoordinates()
        {
        }
        ;

        XYCoordinates(T iX, T iY) :
            x(iX), y(iY)
        {
        }

        XYCoordinates(T iN) :
            x(iN), y(iN)
        {
        }

        //Copy constructor - can be used to perform type converstion 
        template<class OldTypeT>
        XYCoordinates<T>(const XYCoordinates<OldTypeT> & iOldXYCoordinates)
        {
          x = static_cast<T>(iOldXYCoordinates.x);
          y = static_cast<T>(iOldXYCoordinates.y);
        }

        //Equality
        bool operator==(const XYCoordinates<T> right)
        {
          if (x != right.x)
          {
            return false;
          }
          if (y != right.y)
          {
            return false;
          }
          return true;
        }

        //Vector addition
        XYCoordinates<T> operator+(const XYCoordinates<T> right) const
        {
          return XYCoordinates(x + right.x, y + right.y);
        }

        //Vector addition
        XYCoordinates<T>& operator+=(const XYCoordinates<T> right)
        {
          x += right.x;
          y += right.y;

          return *this;
        }

        //Vector subtraction
        XYCoordinates<T> operator-(const XYCoordinates<T> right) const
        {
          return XYCoordinates(x - right.x, y - right.y);
        }

        //Scalar multiplication
        template<class MultiplierT>
        XYCoordinates<T> operator*(const MultiplierT multiplier) const
        {
          return XYCoordinates(x * multiplier, y * multiplier);
        }

        //Updates the XYCoordinates with the smallest of each
        //of the x and y co-ordinatess independently of both XYCoordinates
        void UpdatePointwiseMin(const XYCoordinates<T>& iCompareLocation)
        {
          x = util::NumericalFunctions::min(x, iCompareLocation.x);

          y = util::NumericalFunctions::min(y, iCompareLocation.y);
        }

        //Updates the XYCoordinates with the largest of each
        //of the x and y co-ordinates independently of both XYCoordinates
        void UpdatePointwiseMax(const XYCoordinates<T>& iCompareLocation)
        {
          x = util::NumericalFunctions::max(x, iCompareLocation.x);

          y = util::NumericalFunctions::max(y, iCompareLocation.y);
        }

        static XYCoordinates<T> MaxLimit()
        {
          return XYCoordinates(std::numeric_limits<T>::max());
        }

        static XYCoordinates<T> MinLimit()
        {
          return XYCoordinates(std::numeric_limits<T>::min());
        }
    };

    template<class T>
    XYCoordinates<T> operator+(const XYCoordinates<T> left, const XYCoordinates<T> right)
    {
      return left + right;
    }

    template<class T>
    XYCoordinates<T> operator-(const XYCoordinates<T> left, const XYCoordinates<T> right)
    {
      return left - right;
    }

    template<class TLeft, class TRight>
    XYCoordinates<TRight> operator*(const TLeft left, const XYCoordinates<TRight> right)
    {
      return right * left;
    }
  }
}

#endif // HEMELB_VIS_XYCOORDINATES_H
