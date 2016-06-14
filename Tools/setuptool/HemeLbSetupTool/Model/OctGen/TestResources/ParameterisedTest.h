// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_TESTRESOURCES_PARAMETERISEDTEST_H
#define HEMELBSETUPTOOL_TESTRESOURCES_PARAMETERISEDTEST_H

#include <cppunit/TestCase.h>

template <class FixtureT, class ArgT>
class ParameterisedTest : public CppUnit::TestCase {
public:
  typedef void (FixtureT::*TestMethod)(ArgT);
  ParameterisedTest(std::string name, FixtureT* fix, TestMethod f, ArgT a) :
    CppUnit::TestCase(name), fixture(fix), func(f), arg(a) {
  }
  ParameterisedTest(const ParameterisedTest* other) = delete;
  ParameterisedTest& operator=(const ParameterisedTest& other) = delete;
  
  void runTest() {
    (fixture->*func)(arg);
  }
  void setUp() { 
    fixture->setUp(); 
  }
  void tearDown() { 
    fixture->tearDown(); 
  }
private:
  FixtureT* fixture;
  TestMethod func;
  ArgT arg;
};

#define PARAMETERISED_TEST(Method, ParamT, Param)			\
  CPPUNIT_TEST_SUITE_ADD_TEST((new ParameterisedTest<TestFixtureType, ParamT>(context.getTestNameFor(#Method #Param), \
									     context.makeFixture(), \
									     &TestFixtureType::Method, \
									      Param)))

#endif
