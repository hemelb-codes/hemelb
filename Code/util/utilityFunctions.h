#ifndef HEMELB_UTILITYFUNCTIONS_H
#define HEMELB_UTILITYFUNCTIONS_H

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
    };

    class NumericalMethods
    {
      public:
        /*
         * The Newton_Raphson method takes in a functor with the operator () overloaded, which
         * should return void and take in three doubles: x, f, df. x is the variable to
         * be solved for, f and df are the function and derivative values at x respectively.
         * The function calculates f and df at the given x and stores in the given f and df.
         * The other two arguments for NewtonRaphson are the initial guess and desired accuracy.
         */
        template<class F>
        static double NewtonRaphson(F* func, double x0, double alphaAcc)
        {
          double x = x0, dx;
          double f, df;

          for (int i = 0; i < 20; i++)
          {
            (*func)(x, f, df);

            dx = f / df;
            x -= dx;

            if (fabs(dx) < alphaAcc)
            {
              return x;
            }
          }

          /*
           * TODO: Implement some sensible way of dealing with too many iterations
           */

          return x;
        }

        template<class F>
        static double Brent(F* func, double xl, double xh, double alphaAcc, double fAcc)
        {
          double a = xl, fa;
          double b = xh, fb;
          double c = a, fc;
          double d; // First set after first iteration hence mflag
          double s, fs;

          (*func)(a, fa);
          (*func)(b, fb);
          fc = fa;
          fs = fb;

          // The task of verifying whether a root is enclosed is left to caller

          if (fabs(fa) < fabs(fb))
          {
            double temp = fa;
            fa = fb;
            fb = temp;
            temp = a;
            a = b;
            b = temp;
          }

          bool mflag = true;

          while (fabs(b - a) > alphaAcc && fabs(fb) > fAcc && fabs(fs) > fAcc)
          {
            if (fa != fc && fb != fc)
            {
              s = (a * fb * fc) / ( (fa - fb) * (fa - fc)) + (b * fa * fc) / ( (fb - fa)
                  * (fb - fc)) + (c * fa * fb) / ( (fc - fa) * (fc - fb));
            }
            else
            {
              s = b - fb * (b - a) / (fb - fa);
            }

            if ( (a < b && s < (3 * a + b) / 4.0 && s > b) || (a > b && s > (3 * a + b) / 4.0 && s
                < b) || (mflag && fabs(s - b) >= fabs(b - c) / 2.0) || (!mflag && fabs(s - b)
                >= fabs(c - d) / 2.0) || (mflag && fabs(b - c) < alphaAcc) || (!mflag
                && fabs(c - d) < alphaAcc))
            {
              s = (a + b) / 2.0;
              mflag = true;
            }
            else
            {
              mflag = false;
            }

            (*func)(s, fs);
            d = c;
            c = b;
            fc = fb;

            if (fa * fs < 0)
            {
              b = s;
              fb = fs;
            }
            else
            {
              a = s;
              fa = fs;
            }

            if (fabs(fa) < fabs(fb))
            {
              double temp = fa;
              fa = fb;
              fb = temp;
              temp = a;
              a = b;
              b = temp;
            }
          }

          if (fabs(fb) < fabs(fs))
            return b;
          else
            return s;
        }
    };

    // Returns the number of seconds to 6dp elapsed since the Epoch
    double myClock();
  }
}

#endif // HEMELB_UTILITYFUNCTIONS_H
