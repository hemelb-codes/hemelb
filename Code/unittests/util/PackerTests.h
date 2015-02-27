//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_UTIL_PACKER_H
#define HEMELB_UNITTESTS_UTIL_PACKER_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "util/Packer.h"

namespace hemelb
{
  namespace unittests
  {
    namespace util
    {
      class PackerTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(PackerTests);
          CPPUNIT_TEST(TestInteger);
          CPPUNIT_TEST(TestDouble);
          CPPUNIT_TEST(TestVector3D);
          CPPUNIT_TEST(TestVector);
          CPPUNIT_TEST_SUITE_END();

        public:

          void TestInteger()
          {
            int64_t const i0 = 20000000;
            int16_t const j0 = 2;
            int64_t i = i0;
            int16_t j = j0;
            hemelb::util::Packer send;
            send << i << j;
            CPPUNIT_ASSERT_EQUAL(sizeof(int64_t) + sizeof(int16_t), send.size());
            i = 0; j = 0;

            hemelb::util::Packer receive(send);
            receive.reset_read();
            CPPUNIT_ASSERT(i0 != i);
            CPPUNIT_ASSERT(j0 != j);
            receive >> i >> j;
            CPPUNIT_ASSERT_EQUAL(i0, i);
            CPPUNIT_ASSERT_EQUAL(j0, j);
          }
          void TestDouble()
          {
            double const d0 = 3.14;
            float const f0 = 3.1;
            double d = d0;
            float f = f0;
            hemelb::util::Packer send;
            send << d << f;
            send.reset_read();
            d = 0e0; f = 0e0;
            send >> d >> f;
            CPPUNIT_ASSERT_DOUBLES_EQUAL(d0, d, 1e-16);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(f0, f, 1e-8);
          }
          void TestVector3D()
          {
            hemelb::util::Vector3D<int> const v0 = {0, 1, 2};
            hemelb::util::Vector3D<int> v = {2, 1, 0};
            hemelb::util::Packer send;
            send << v0;
            send.reset_read();
            send >> v;
            CPPUNIT_ASSERT_EQUAL(v0.x, v.x);
            CPPUNIT_ASSERT_EQUAL(v0.y, v.y);
            CPPUNIT_ASSERT_EQUAL(v0.z, v.z);
          }
          void TestVector()
          {
            std::vector<hemelb::util::Vector3D<int>> const v0 = {{0, 1, 2}, {2, 3, 4}, {5, 6, 7}};
            std::vector<hemelb::util::Vector3D<int>> v = v0;
            hemelb::util::Packer send;
            send << v0;
            send.reset_read();
            send >> v;
            for(int i = 0; i < v0.size(); ++i)
            {
              CPPUNIT_ASSERT_EQUAL(v0[i].x, v[i].x);
              CPPUNIT_ASSERT_EQUAL(v0[i].y, v[i].y);
              CPPUNIT_ASSERT_EQUAL(v0[i].z, v[i].z);
            }
          }

          void TestMap()
          {
            std::map<int, std::string> const m0 = {{0, "zero"}, {1, "one"}};
            std::map<int, std::string> m;
            hemelb::util::Packer packer;
            packer << m0;
            packer.reset_read();
            packer >> m;
            CPPUNIT_ASSERT_EQUAL(m0.size(), m.size());
            for(auto const &element: m0)
            {
              CPPUNIT_ASSERT_EQUAL(1ul, m.count(element.first));
              CPPUNIT_ASSERT_EQUAL(element.second, m[element.first]);
            }
          }
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(PackerTests);

    }
  }
}

#endif /* HEMELB_UNITTESTS_UTIL_MATRIX3DTESTS_H */
