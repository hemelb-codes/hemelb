#ifndef TESTS_H
#define TESTS_H

#include "lb/GlobalLatticeData.h"

class Tests
{
  public:
    static void RunTests();

  private:
    static void TestDomainDecomposition();

    static void MakeAGlobLatDat(hemelb::lb::GlobalLatticeData & oGlobLatDat);
};

#endif /* TESTS_H */
