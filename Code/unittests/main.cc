// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#define HEMELB_DOING_UNITTESTS
#include <cppunit/XmlOutputter.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/extensions/HelperMacros.h>
#include <stdexcept>
#include <unistd.h>

#include "unittests/helpers/helpers.h"
#include "debug/Debugger.h"
#include HEMELB_UNITTEST_INCLUDE

int main(int argc, char **argv)
{
  // Start MPI and the logger.
  hemelb::net::MpiEnvironment mpi(argc, argv);
  hemelb::log::Logger::Init();

  hemelb::net::MpiCommunicator commWorld = hemelb::net::MpiCommunicator::World();

  // Read options
  std::ostream * reportto = &std::cerr;
  std::ofstream reportfile;
  int opt;
  bool debug = false;

  while ( (opt = getopt(argc, argv, "o:d")) != -1)
  {
    switch (opt)
    {
      case 'o':
        reportfile.open(optarg);
        reportto = &reportfile;
        break;
      case 'd':
        debug = true;
        break;
    }
  }
  // Start the debugger (no-op if HEMELB_USE_DEBUGGER is OFF)
  hemelb::debug::Debugger::Init(debug, argv[0], commWorld);

  // Initialise the global IOCommunicator.
  hemelb::net::IOCommunicator testCommunicator(commWorld);
  hemelb::unittests::helpers::HasCommsTestFixture::Init(testCommunicator);

  std::string testPath = (optind < argc) ?
    std::string(argv[optind]) :
    "";
  // Create the event manager and test controller
  CppUnit::TestResult controller;

  // Add a listener that collects test result
  CppUnit::TestResultCollector result;
  controller.addListener(&result);

  // Add a listener that print dots to stdout as test run.
  CppUnit::BriefTestProgressListener progress;
  controller.addListener(&progress);

  // Add the top suite to the test runner
  CppUnit::TestRunner runner;

  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  runner.addTest(registry.makeTest());

  try
  {
    std::cout << "Running " << testPath;
    runner.run(controller, testPath);
    // Print test XML output to stderr
    CppUnit::XmlOutputter outputter(&result, *reportto);
    outputter.write();
  }
  catch (std::invalid_argument &e) // Test path not resolved
  {
    std::cerr << std::endl << "ERROR: " << e.what() << std::endl;
    reportfile.close();
    return 1;
  }

  reportfile.close();
  return result.wasSuccessful() ?
    0 :
    1;
}

