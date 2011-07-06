#ifndef HEMELB_UTILITYFUNCTIONS_H
#define HEMELB_UTILITYFUNCTIONS_H

#include "util/utilityFunctions.h"
#include <math.h>
#include <iostream>
#include <cstdlib>

// Static class for simple functions that could be useful in many places
namespace hemelb
{
  namespace util
  {
    class NumericalFunctions
    {
      public:
        // Return the smaller of the two numbers
        template<typename T>
        static T min(T a, T b)
        {
          if (a < b)
          {
            return a;
          }
          else
          {
            return b;
          }
        }

        //Return the larger of the two numbers
        template<typename T>
        static T max(T a, T b)
        {
          if (a > b)
          {
            return a;
          }
          else
          {
            return b;
          }
        }

        // If number < lowerBound, returns lowerBound. If number >
        // upperBound, returns upperBound If number between bounds, returns
        // number.  Consider the behaviour undefined if lowerBound >
        // upperBound, so don't try it!
        template<typename T>
        static T enforceBounds(T number, T lowerBound, T upperBound)
        {
          return max<T> (lowerBound, min(number, upperBound));
        }

        /*
         * The Newton_Raphson method takes in a functor with the operator () overloaded, which
         * should return void and take in three doubles: x, f, df. x is the variable to
         * be solved for, f and df are the function and derivative values at x respectively.
         * The function calculates f and df at the given x and stores in the given f and df.
         * The other two arguments for NewtonRaphson are the initial guess and desired accuracy.
         */
        template<class F>
        static double NewtonRaphson(F* func, double x0, double acc)
        {
          double x = x0;
          double f, df;

          for (int i = 0; i < 20; i++)
          {
            (*func)(x, f, df);
            x0 = x;
            x = x0 - (f / df);
            if (fabs(x - x0) < acc)
            {
              return x;
            }
          }

          /*
           * TODO: Implement some sensible way of dealing with too many iterations
           */

          return x;
        }
    };

    // Returns the number of seconds to 6dp elapsed since the Epoch
    double myClock();
  }
}

#endif // HEMELB_UTILITYFUNCTIONS_H
