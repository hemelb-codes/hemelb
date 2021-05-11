// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_LOADINGTIMMMESHTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_LOADINGTIMMMESHTESTS_H

#include <sstream>
#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "redblood/Mesh.h"
#include "redblood/MeshIO.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      // Tests we can load .msh files output by Timm's code
      class LoadingTimmMeshTests : public BasisFixture
      {
          CPPUNIT_TEST_SUITE (LoadingTimmMeshTests);
          CPPUNIT_TEST (testWrongIcoFile);
          CPPUNIT_TEST (testCorrectedIcoFile);
          CPPUNIT_TEST (testSimulationOutputFile);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
          }

          bool checkMeshOrientation(MeshData const &mesh)
          {
            MeshData copy(mesh);
            orientFacets(copy);
            for (size_t i(0); i < mesh.facets.size(); ++i)
            {
              if (mesh.facets[i] != copy.facets[i])
              {
                return false;
              }
            }
            return true;
          }

          std::shared_ptr<MeshData> readMeshFromFileNoFix(const std::string filename)
          {
            auto fullPath = resources::Resource(filename).Path();
            CPPUNIT_ASSERT(util::file_exists(fullPath.c_str()));
            return redblood::KruegerMeshIO{}.readFile(fullPath, false);
          }

          void testWrongIcoFile()
          {
            // This is one of Timm's original ico files and has a face with the wrong vertex orientation leading to an inward pointing normal
            // Loading it without fixing face orientation should lead to inconsistent orientation
            auto mesh = readMeshFromFileNoFix("rbc_ico_720.msh");
            CPPUNIT_ASSERT(!checkMeshOrientation(*mesh));
          }

          void testCorrectedIcoFile()
          {
            // This is the same ico file with the face fixed
            // Loading it without fixing face orientation should be fine
            auto mesh = readMeshFromFileNoFix("rbc_ico_720_correct.msh");
            CPPUNIT_ASSERT(checkMeshOrientation(*mesh));
          }

          void testSimulationOutputFile()
          {
            // This is a mesh taken from a simulation in Timm's code, which should have correct face orientation by definition
            auto mesh = readMeshFromFileNoFix("992Particles_rank3_26_t992.msh");
            CPPUNIT_ASSERT(checkMeshOrientation(*mesh));
          }
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (LoadingTimmMeshTests);
    }
  }
}

#endif  // ONCE
