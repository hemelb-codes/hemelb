// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_MPI_H
#define HEMELB_UNITTESTS_REDBLOOD_MPI_H
#define HEMELB_DOING_UNITTESTS
#include <unittests/redblood/parallel/CellParallelizationTests.h>
#include <unittests/redblood/parallel/MPISpreadForcesTests.h>
#include <unittests/redblood/parallel/MPIIntegrateVelocities.h>
#include <unittests/redblood/parallel/ParallelFixtureTests.h>
#include <unittests/redblood/parallel/MPILockStepTests.h>
#include <unittests/redblood/parallel/GraphCommsTests.h>
#include <unittests/redblood/parallel/RBCInserterWithPerturbationParallelTests.h>
#undef HEMELB_DOING_UNITTESTS
#endif
