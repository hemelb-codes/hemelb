#include "timing.h"


double myClock ()
{
  return (double)clock () * (1.0 / (double)CLOCKS_PER_SEC);
}
