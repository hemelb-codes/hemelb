// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <cppunit/XmlOutputter.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/extensions/HelperMacros.h>
#include <stdexcept>
#include "unittests/helpers/helpers.h"
#include "unittests/lbtests/lbtests.h"
#include "unittests/vistests/vistests.h"
#include "unittests/io/io.h"
#include "unittests/reporting/reporting.h"
#include "unittests/configuration/configuration.h"
#include "unittests/geometry/geometry.h"
#include "unittests/SimulationMasterTests.h"
#include "unittests/extraction/extraction.h"
#include "unittests/net/net.h"
#include "unittests/multiscale/multiscale.h"
#ifdef HEMELB_BUILD_MULTISCALE
  #include "unittests/multiscale/mpwide/mpwide.h"
#endif
#include "unittests/util/util.h"
#include <unistd.h>

int main(int argc, char **argv)
{
  hemelb::net::MpiEnvironment mpi(argc, argv);
  hemelb::log::Logger::Init();

  hemelb::net::MpiCommunicator testCommunicator = hemelb::net::MpiCommunicator::World();
  hemelb::net::NetworkTopology::Instance()->Init(testCommunicator);

  std::ostream * reportto=&std::cerr;
  std::ofstream reportfile;
  int opt;
  while((opt=getopt(argc,argv,"o:"))!=-1){
    switch (opt) {
      case 'o':
        reportfile.open(optarg);
        reportto=&reportfile;
        break;
    }
  }

  std::string testPath = (optind < argc)
    ? std::string(argv[optind])
    : "";
  // Create the event manager and test controller
  CppUnit::TestResult controller;

  // Add a listener that colllects test result
  CppUnit::TestResultCollector result;
  controller.addListener(&result);

  // Add a listener that print dots to stdout as test run.
  CppUnit::BriefTestProgressListener progress;
  controller.addListener(&progress);

  // Add the top suite to the test runner
  CppUnit::TestRunner runner;

  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  runner.addTest( registry.makeTest() );

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
  return result.wasSuccessful()
    ? 0
    : 1;
}

