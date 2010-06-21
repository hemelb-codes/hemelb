// Static class for simple functions that could be useful in many places
class UtilityFunctions
{
  public:
    // Simple integer comparisons.
    static int min (int a, int b);
    static int max (int a, int b);
    static int enforceBounds(int number, int lowerBound, int upperBound);

    // Returns the number of seconds to 6dp elapsed since the Epoch
    static double myClock ();
};