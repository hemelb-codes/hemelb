#ifndef HEMELB_UTILITYFUNCTIONS_H
#define HEMELB_UTILITYFUNCTIONS_H

#include "util/utilityFunctions.h"
#include "log/Logger.h"
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <vector>

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

        template<typename T>
        static T LinearInterpolate(std::vector<T> &xs, std::vector<T> &ys, T x)
        {
          int i = 0;

          if (hemelb::log::Logger::ShouldDisplay<hemelb::log::Debug>())
          {
            if (x < xs[0] || x > xs[xs.size() - 1])
            {
              log::Logger::Log<log::Info, log::OnePerCore>("Linear Interpolation beyond bounds: ");
              log::Logger::Log<log::Info, log::OnePerCore>("%f is not between %f and %f",
                                                           x,
                                                           xs[0],
                                                           xs[xs.size() - 1]);
            }
          }

          while (! (x >= xs[i] && x <= xs[i + 1]))
          {
            i++;
          }

          // In case repeat values are read in or discontinuities are present in the trace
          // Average taken in case of a discontinuity. Repeat points should just contract to one.
          //if (hemelb::log::Logger::ShouldDisplay<hemelb::log::Debug>())
          {
            if (xs[i] == xs[i + 1])
            {
              log::Logger::Log<log::Info, log::OnePerCore>("Multiple points for same x value in LinearInterpolate: ");
              log::Logger::Log<log::Info, log::OnePerCore>("(%f, %f) and (%f, %f). Division by zero!",
                                                           xs[i],
                                                           ys[i],
                                                           xs[i + 1],
                                                           ys[i + 1]);
            }
          }

          // Linear interpolation of function f(x) between two points A and B
          // f(A) + (fraction along x axis between A and B) * (f(B) - f(A))
          return (ys[i] + (x - xs[i]) / (xs[i + 1] - xs[i]) * (ys[i + 1] - ys[i]));
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

        /*
         * Finds root using Brent's method. Needs to be given a bracket enclosing the root.
         * The caller must check if a root is enclosed so that he can specify the result in that case
         * Since it must check for this it will have the values of the function at those points
         * so they need to be passed on as well.
         */
        template<class F>
        static double Brent(F* func,
                            double xl,
                            double fl,
                            double xh,
                            double fh,
                            double alphaAcc,
                            double fAcc)
        {
          double a = xl, fa = fl;
          double b = xh, fb = fh;
          double c = a, fc = fa;
          double d; // First set after first iteration hence mflag
          double s, fs = fb;

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
