
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_EXTRACTION_LOCALPROPERTYOUTPUTTESTS_H
#define HEMELB_UNITTESTS_EXTRACTION_LOCALPROPERTYOUTPUTTESTS_H

#include <string>
#include <cstdio>

#include <cppunit/TestFixture.h>

#include "io/formats/extraction.h"
#include "io/writers/xdr/XdrMemReader.h"
#include "extraction/PropertyOutputFile.h"
#include "extraction/OutputField.h"
#include "extraction/WholeGeometrySelector.h"

#include "unittests/helpers/HasCommsTestFixture.h"
#include "unittests/extraction/DummyDataSource.h"

namespace hemelb
{
  namespace unittests
  {
    namespace extraction
    {
      class LocalPropertyOutputTests : public helpers::HasCommsTestFixture
      {
          CPPUNIT_TEST_SUITE (LocalPropertyOutputTests);
          CPPUNIT_TEST (TestStringWrittenLength);
          CPPUNIT_TEST (TestWrite);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            helpers::HasCommsTestFixture::setUp();
            epsilon = 1e-5;

            simpleOutFile.filename = tempOutFileName;
            // The code won't overwrite any existing file
            std::remove(tempOutFileName);

            simpleOutFile.frequency = 100;

            simpleOutFile.geometry = new hemelb::extraction::WholeGeometrySelector();

            hemelb::extraction::OutputField pressure;
            pressure.name = "Pressure";
            pressure.type = hemelb::extraction::OutputField::Pressure;
            simpleOutFile.fields.push_back(pressure);

            hemelb::extraction::OutputField velocity;
            velocity.name = "Velocity";
            velocity.type = hemelb::extraction::OutputField::Velocity;
            simpleOutFile.fields.push_back(velocity);

            simpleDataSource = new DummyDataSource();

            writtenFile = NULL;
            propertyWriter = NULL;
            writtenMainHeader = new char[hemelb::io::formats::extraction::MainHeaderLength];

            fieldHeaderLength = 0x30;
            writtenFieldHeader = new char[fieldHeaderLength];

          }
          void tearDown()
          {
            delete simpleDataSource;

            delete[] writtenMainHeader;
            delete[] writtenFieldHeader;

            // Delete the LocalPropertyOutput
            if (propertyWriter != NULL)
              delete propertyWriter;
            // Close the file and delete it
            if (writtenFile != NULL)
              std::fclose(writtenFile);
            std::remove(tempOutFileName);

            helpers::HasCommsTestFixture::tearDown();
          }

          void TestStringWrittenLength()
          {
            // Check the zero length string
            std::string s0;
            CPPUNIT_ASSERT_EQUAL(size_t(4), hemelb::io::formats::extraction::GetStoredLengthOfString(s0));

            // This should have no padding
            std::string s1("Fish");
            CPPUNIT_ASSERT_EQUAL(size_t(8), hemelb::io::formats::extraction::GetStoredLengthOfString(s1));

            // This must be padded up to 8 bytes
            std::string s2("A");
            CPPUNIT_ASSERT_EQUAL(size_t(8), hemelb::io::formats::extraction::GetStoredLengthOfString(s2));

          }

          void TestWrite()
          {
            // Create the writer object; this should write the headers.
            propertyWriter = new hemelb::extraction::LocalPropertyOutput(*simpleDataSource, &simpleOutFile, Comms());

            // Open the file
            writtenFile = std::fopen(simpleOutFile.filename.c_str(), "r");

            // Assert that the file is there
            CPPUNIT_ASSERT(writtenFile != NULL);

            // Read the main header.
            size_t nRead = std::fread(writtenMainHeader,
                                      1,
                                      hemelb::io::formats::extraction::MainHeaderLength,
                                      writtenFile);
            // Check we read enough
            CPPUNIT_ASSERT_EQUAL(size_t(hemelb::io::formats::extraction::MainHeaderLength), nRead);

            // Check it matches
            const char expectedMainHeader[] = "\x68\x6C\x62\x21"
                "\x78\x74\x72\x04"
                "\x00\x00\x00\x04"
                "\x3F\x33\xA9\x2A"
                "\x30\x55\x32\x61"
                "\x3F\xA1\x68\x72"
                "\xB0\x20\xC4\x9C"
                "\x3F\x50\x62\x4D"
                "\xD2\xF1\xA9\xFC"
                "\x3F\xB2\xF1\xA9"
                "\xFB\xE7\x6C\x8B"
                "\x00\x00\x00\x00"
                "\x00\x00\x00\x40"
                "\x00\x00\x00\x02"
                "\x00\x00\x00\x30";
            for (int i = 0; i < hemelb::io::formats::extraction::MainHeaderLength; ++i)
            {
              CPPUNIT_ASSERT_EQUAL(expectedMainHeader[i], writtenMainHeader[i]);
            }

            // Read the field header.
            nRead = std::fread(writtenFieldHeader, 1, fieldHeaderLength, writtenFile);

            // Check we read enough
            CPPUNIT_ASSERT_EQUAL(fieldHeaderLength, nRead);

            const char expectedFieldHeader[] = "\x00\x00\x00\x08"
                "\x50\x72\x65\x73"
                "\x73\x75\x72\x65"
                "\x00\x00\x00\x01"
                "\x00\x00\x00\x00"
                "\x00\x00\x00\x00"
                "\x00\x00\x00\x08"
                "\x56\x65\x6C\x6F"
                "\x63\x69\x74\x79"
                "\x00\x00\x00\x03"
                "\x00\x00\x00\x00"
                "\x00\x00\x00\x00";

            for (size_t i = 0; i < fieldHeaderLength; ++i)
            {
              CPPUNIT_ASSERT_EQUAL(expectedFieldHeader[i], writtenFieldHeader[i]);
            }

            // Now going to write the body.
            // Create some pseudo-random data
            simpleDataSource->FillFields();
            // Write it
            propertyWriter->Write(0);

            CheckDataWriting(simpleDataSource, 0, writtenFile);

            // Get some new data
            simpleDataSource->FillFields();
            // This should NOT write
            propertyWriter->Write(10);
            // This SHOULD write
            propertyWriter->Write(100);

            // The previous call to CheckDataWriting() sets the EOF indicator in writtenFile,
            // the previous call to Write() ought to unset it but it isn't working properly in
            // my current version of libc.
            std::clearerr(writtenFile);

            CheckDataWriting(simpleDataSource, 100, writtenFile);
          }

        private:
          void CheckDataWriting(DummyDataSource* datasource, uint64_t timestep, FILE* file)
          {
            // The file should have an entry for each lattice point, consisting
            // of 3D grid coords, pressure (with an offset of 80) and 3D velocity.
            // This gives 3*4 + 4 + 3*4 = 28 bytes per site.
            long siteCount = 0;
            datasource->Reset();
            while (datasource->ReadNext())
            {
              ++siteCount;
            }

            // We also have the iteration number, a long
            size_t expectedSize = 8 + 28 * siteCount;

            // Attempt to read one extra byte, to make sure we aren't under-reading
            char* contentsBuffer = new char[expectedSize];
            size_t nRead = std::fread(contentsBuffer, 1, expectedSize + 1, file);

            CPPUNIT_ASSERT_EQUAL(expectedSize, nRead);

            // Create an XDR buffer reader to process the file for comparison
            hemelb::io::writers::xdr::XdrMemReader reader(contentsBuffer, expectedSize);

            // The timestep should be correct
            uint64_t readTimestep;
            reader.readUnsignedLong(readTimestep);

            CPPUNIT_ASSERT_EQUAL(timestep, readTimestep);

            // Now iterate over the sites.
            datasource->Reset();
            while (datasource->ReadNext())
            {
              // Read the grid, which should be the same
              LatticeVector grid = datasource->GetPosition();
              unsigned x, y, z;
              reader.readUnsignedInt(x);
              reader.readUnsignedInt(y);
              reader.readUnsignedInt(z);

              CPPUNIT_ASSERT_EQUAL((unsigned) grid.x, x);
              CPPUNIT_ASSERT_EQUAL((unsigned) grid.y, y);
              CPPUNIT_ASSERT_EQUAL((unsigned) grid.z, z);

              // Read the pressure, which should be an offset of the
              // reference pressure away.
              float pressure;
              reader.readFloat(pressure);

              CPPUNIT_ASSERT_DOUBLES_EQUAL(datasource->GetPressure(), (REFERENCE_PRESSURE_mmHg + (double) pressure), epsilon);

              // Read the velocity and compare
              PhysicalVelocity velocity = datasource->GetVelocity();
              float vx, vy, vz;
              reader.readFloat(vx);
              reader.readFloat(vy);
              reader.readFloat(vz);

              CPPUNIT_ASSERT_DOUBLES_EQUAL(velocity.x, (double) vx, epsilon);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(velocity.y, (double) vy, epsilon);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(velocity.z, (double) vz, epsilon);
            }

            delete[] contentsBuffer;
          }

          hemelb::extraction::PropertyOutputFile simpleOutFile;
          DummyDataSource* simpleDataSource;

          hemelb::extraction::LocalPropertyOutput* propertyWriter;
          double epsilon;
          FILE* writtenFile;
          char* writtenMainHeader;
          size_t fieldHeaderLength;
          char* writtenFieldHeader;
          static const char* tempOutFileName;
      };
      const char* LocalPropertyOutputTests::tempOutFileName = "simple.xtr";
      CPPUNIT_TEST_SUITE_REGISTRATION (LocalPropertyOutputTests);
    }
  }
}

#endif // HEMELB_UNITTESTS_EXTRACTION_LOCALPROPERTYOUTPUTTESTS_H
