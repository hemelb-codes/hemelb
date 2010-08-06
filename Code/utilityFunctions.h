// Static class for simple functions that could be useful in many places
namespace util
{
  // Simple integer comparisons.
  int min (int a, int b);
  int max (int a, int b);
  int enforceBounds(int number, int lowerBound, int upperBound);
  
  // Returns the number of seconds to 6dp elapsed since the Epoch
   double myClock ();
};
