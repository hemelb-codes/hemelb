#include "timing.h"


double myClock ()
{
  return (double)clock () * (1. / (double)CLOCKS_PER_SEC);
}
