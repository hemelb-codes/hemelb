// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

// In hemelb/Code/tests/extraction/LocalPropertyHdf5OutputTests.cc

#include <catch2/catch.hpp>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>

// HemeLB includes
#include "extraction/LocalPropertyHdf5Output.h"
#include "extraction/PropertyOutputFile.h"
#include "extraction/OutputField.h"
#include "extraction/WholeGeometrySelector.h"
#include "tests/helpers/HasCommsTestFixture.h"
#include "tests/extraction/DummyDataSource.h"
#include "net/IOCommunicator.h"

// HDF5 C API include
#include <hdf5.h>

namespace hemelb
{
  namespace tests
  {
    void CheckHdf5Output(const std::vector<double>& expected_pressures, 
                         const extraction::PropertyOutputFile& spec, 
                         const net::IOCommunicator& comms, 
                         unsigned long timestep)
    {
      if (comms.Rank() != 0) return;
      
      std::string h5_filename = spec.filename.string();

      hid_t file_id = H5Fopen(h5_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      INFO("Checking if HDF5 file was created and can be opened: " << h5_filename);
      REQUIRE(file_id >= 0);

      std::string group_name = "step_" + std::to_string(timestep);
      hid_t group_id = H5Gopen(file_id, group_name.c_str(), H5P_DEFAULT);
      INFO("Checking if group exists: " << group_name);
      REQUIRE(group_id >= 0);

      hid_t dset_id = H5Dopen(group_id, "pressure", H5P_DEFAULT);
      INFO("Checking if dataset 'pressure' exists in group.");
      REQUIRE(dset_id >= 0);
      
      std::vector<double> read_pressures(expected_pressures.size());
      H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_pressures.data());
      
      REQUIRE(read_pressures.size() == expected_pressures.size());
      
      bool mismatch_found = false;
      for (size_t i = 0; i < read_pressures.size(); ++i) {
          if (read_pressures[i] != Approx(expected_pressures[i])) {
              if (!mismatch_found) { // only print mismatch details once
                  std::cout << "\n--- Mismatch Details for 'pressure' ---\n";
                  mismatch_found = true;
              }
              std::cout << "[HDF5] read[" << i << "] = " << read_pressures[i]
                        << ", expected = " << expected_pressures[i] << std::endl;
          }
      }
      // if mismatch was found, print a summary
      REQUIRE_THAT(read_pressures, Catch::Matchers::Approx(expected_pressures));

      H5Dclose(dset_id);
      H5Gclose(group_id);
      H5Fclose(file_id);
    }


    TEST_CASE_METHOD(helpers::HasCommsTestFixture, "HDF5 Property Output functionality", "[extraction][hdf5]")
    {
      const char* tempH5FileName = "final_test_output.h5";
      const char* tempXmfFileName = "final_test_output.xmf";
      
      if (Comms().Rank() == 0) {
        std::remove(tempH5FileName);
        std::remove(tempXmfFileName);
      }
      Comms().Barrier();
      
      // use DummyDataSource for testing
      auto dataSource = std::make_unique<DummyDataSource>(42); 
      
      extraction::PropertyOutputFile spec;
      spec.filename = tempH5FileName;
      spec.frequency = 100; // only write every 100 steps
      spec.geometry = util::make_clone_ptr<extraction::WholeGeometrySelector>();

      extraction::OutputField pressure_field;
      pressure_field.name = "pressure";
      pressure_field.src = extraction::source::Pressure{};
      spec.fields.push_back(pressure_field);
      
      std::vector<double> ground_truth_pressures;

      {
        auto hdf5Writer = std::make_unique<extraction::LocalPropertyHdf5Output>(*dataSource, spec, Comms());
        
        for (unsigned long step = 0; step <= 100; ++step) {
            if (hdf5Writer->ShouldWrite(step)) {
                dataSource->FillFields();
                if (step == 100) {
                    ground_truth_pressures.clear();
                    dataSource->Reset();
                    while(dataSource->ReadNext()){
                        ground_truth_pressures.push_back(dataSource->GetPressure());
                    }
                }
            }
            hdf5Writer->Write(step, 100);
        }
      }

      REQUIRE(!ground_truth_pressures.empty());
      CheckHdf5Output(ground_truth_pressures, spec, Comms(), 100);

      if (Comms().Rank() == 0) {
        std::remove(tempH5FileName);
        std::remove(tempXmfFileName);
      }
    }
  }
}