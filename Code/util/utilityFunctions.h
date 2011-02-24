#ifndef HEMELB_UTILITYFUNCTIONS_H
#define HEMELB_UTILITYFUNCTIONS_H

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

    // Returns the number of seconds to 6dp elapsed since the Epoch
    double myClock();
  }
}

#endif // HEMELB_UTILITYFUNCTIONS_H
