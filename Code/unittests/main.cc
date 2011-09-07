#include <cppunit/ui/text/TestRunner.h>

#include "unittests/lbtests/lbtests.h"

int main(int argc, char **argv)
{
  CppUnit::TextUi::TestRunner runner;

  runner.addTest(new hemelb::unittests::lbtests::LbTestSuite());

  runner.run();

  return 0;
}
