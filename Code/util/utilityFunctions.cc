#include "mpiInclude.h"
#include "util/utilityFunctions.h"
#include <math.h>

namespace hemelb
{
  namespace util
  {
    // Returns the number of seconds to 6dp elapsed since the Epoch
    double myClock()
    {
      return MPI_Wtime();
    }

    /*
     * The Newton_Raphson method takes in a functor with the operator () overloaded, which
     * should return void and take in three doubles: x, f, df. x is the variable to
     * be solved for, f and df are the function and derivative values at x respectively.
     * The function calculates f and df at the given x and stores in the given f and df.
     * The other two arguments for NewtonRaphson are the initial guess and desired accuracy.
     */
    template <class F>
    double NewtonRaphson(F &func, double x0, double acc)
    {
      double x = x0;
      double f, df;

      for (int i=0; i < 20; i++)
      {
        func(x, f, df);
        x0 = x;
        x = x0 - (f / df);
        if (fabs(x - x0) < acc)
          return x;
      }

      std::cerr << "Newton Raphson did not converge within allowed number of iterations (20)\n";
      return x;
    }

  }

}
