#ifndef HEMELB_MPIINCLUDE_H
#define HEMELB_MPIINCLUDE_H

#ifndef NOMPI
#ifdef XT3
#include <mpi.h>
#else
#include "mpi.h"
#endif
#endif

#endif // HEMELB_MPIINCLUDE_H
