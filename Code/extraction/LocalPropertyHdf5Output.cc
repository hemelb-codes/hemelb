// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

// LocalPropertyHdf5Output.cc
// In extraction/LocalPropertyHdf5Output.cc
// check hdf5 file, group, separate datasets, and XDMF with mask

#ifdef USE_HDF5

#include "extraction/LocalPropertyHdf5Output.h"
#include "net/IOCommunicator.h"
#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <filesystem>

// check HDF5 error handling
#define H5_CHECK(err) (h5_check((err), __FILE__, __LINE__))
inline herr_t h5_check(herr_t err, const char* file, int line) {
    if (err < 0) {
        fprintf(stderr, "HDF5 Error at %s:%d\n", file, line);
        H5Eprint(H5E_DEFAULT, stderr);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    return err;
}

namespace {
    // parallel write dataset helper function
    void WriteDataset(hid_t loc_id, const std::string& name, const std::vector<double>& buffer, 
                      hsize_t global_rows, hsize_t local_rows, unsigned field_len, 
                      const hemelb::net::IOCommunicator& comms, MPI_Comm mpi_comm) {
        
        hsize_t global_dims[2] = {global_rows, field_len};
        hid_t filespace_id = H5Screate_simple(2, global_dims, NULL);
        H5_CHECK(filespace_id);
        
        hid_t dataset_id = H5Dcreate2(loc_id, name.c_str(), H5T_NATIVE_DOUBLE, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5_CHECK(dataset_id);

        hsize_t local_dims[2] = {local_rows, field_len};
        hid_t memspace_id = H5Screate_simple(2, local_dims, NULL);
        H5_CHECK(memspace_id);

        hsize_t start[2] = {(hsize_t)comms.Scan(local_rows, MPI_SUM) - local_rows, 0};
        hsize_t count[2] = {local_rows, field_len};
        H5_CHECK(H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, start, NULL, count, NULL));
        
        hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
        H5_CHECK(dxpl);
        H5_CHECK(H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE));
        H5_CHECK(H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, dxpl, buffer.data()));
        
        H5_CHECK(H5Pclose(dxpl));
        H5_CHECK(H5Sclose(memspace_id));
        H5_CHECK(H5Dclose(dataset_id));
        H5_CHECK(H5Sclose(filespace_id));
    }
}

namespace hemelb::extraction
{
  LocalPropertyHdf5Output::LocalPropertyHdf5Output(IterableDataSource& dataSource,
                                                     const PropertyOutputFile& outputSpec,
                                                     const net::IOCommunicator& comms)
    : LocalPropertyOutput(dataSource, outputSpec, comms), mpi_comm(comms) {}

  LocalPropertyHdf5Output::~LocalPropertyHdf5Output()
  {
    // before closing the file, write xmf
    if (file_id >= 0) {
      WriteXDMFFile();
      H5Fclose(file_id);
    }
  }
  
  void LocalPropertyHdf5Output::WriteXDMFFile()
  {
    if (comms.Rank() != 0) return;
    
    std::string h5_filename = outputSpec.filename.string();
    std::string xmf_filename = h5_filename;
    size_t pos = xmf_filename.rfind(".h5");
    if (pos != std::string::npos) {
        xmf_filename.replace(pos, 3, ".xmf");
    } else {
        xmf_filename += ".xmf";
    }

    std::filesystem::path h5_path(h5_filename);
    std::string h5_basename = h5_path.filename().string();
    
    std::ofstream xmf_file(xmf_filename);
    xmf_file << "<?xml version=\"1.0\" ?>\n";
    xmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    xmf_file << "<Xdmf Version=\"3.0\">\n";
    xmf_file << "  <Domain>\n";
    xmf_file << "    <Grid Name=\"HemeLB Collection\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

    // Traverse all written timesteps, write each as a grid
    for (unsigned long ts : written_timesteps) {
        xmf_file << "      <Grid Name=\"step_" << ts << "\" GridType=\"Uniform\">\n";
        xmf_file << "        <Time Value=\"" << ts << "\"/>\n";
        xmf_file << "        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << global_site_count << "\"/>\n";
        xmf_file << "        <Geometry GeometryType=\"XYZ\">\n";
        xmf_file << "          <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " 3\" DataType=\"Float\" Precision=\"8\">" << h5_basename << ":/step_" << ts << "/geometry</DataItem>\n";
        xmf_file << "        </Geometry>\n";

        for (const auto& field : outputSpec.fields) {
            unsigned field_len = GetFieldLength(field.src);
            std::string type = (field_len > 1) ? "Vector" : "Scalar";
            xmf_file << "        <Attribute Name=\"" << field.name << "\" AttributeType=\"" << type << "\" Center=\"Node\">\n";
            xmf_file << "          <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " " << field_len << "\" DataType=\"Float\" Precision=\"8\">" << h5_basename << ":/step_" << ts << "/" << field.name << "</DataItem>\n";
            xmf_file << "        </Attribute>\n";
        }
        xmf_file << "      </Grid>\n";
    }

    xmf_file << "    </Grid>\n";
    xmf_file << "  </Domain>\n";
    xmf_file << "</Xdmf>\n";
    xmf_file.close();
  }

  void LocalPropertyHdf5Output::Write(unsigned long timestepNumber, unsigned long totalSteps)
  {
    if (!ShouldWrite(timestepNumber)) return;
    
    // if first time writing, create file
    if (file_id < 0) {
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        H5_CHECK(H5Pset_fapl_mpio(fapl, mpi_comm, MPI_INFO_NULL));
        file_id = H5Fcreate(outputSpec.filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        H5_CHECK(H5Pclose(fapl));
    }
    
    written_timesteps.push_back(timestepNumber);

    // create group for this timestep
    std::string group_name = "step_" + std::to_string(timestepNumber);
    hid_t group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5_CHECK(group_id);
    
    // Get dataset
    std::vector<double> coords_buffer;
    coords_buffer.reserve(local_site_count * 3);
    std::vector<std::vector<double>> field_buffers(outputSpec.fields.size());
    for(size_t i = 0; i < outputSpec.fields.size(); ++i) {
        field_buffers[i].reserve(local_site_count * GetFieldLength(outputSpec.fields[i].src));
    }
    // One traverse get all datasets
    dataSource.Reset();
    while (dataSource.ReadNext()) {
      if (outputSpec.geometry->Include(dataSource, dataSource.GetPosition())) {
        // Get position
        const auto& pos = dataSource.GetPosition();
        coords_buffer.push_back(static_cast<double>(pos.x()));
        coords_buffer.push_back(static_cast<double>(pos.y()));
        coords_buffer.push_back(static_cast<double>(pos.z()));
        // push origin data into field_buffers
        for (size_t i = 0; i < outputSpec.fields.size(); ++i) {
          std::visit([&](auto&& src_type) {
             if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::Pressure>) {
                double val = dataSource.GetPressure();
                field_buffers[i].push_back(val);
                // field_buffers[i].push_back(std::isinf(val) ? -1.0 : val);
             } else if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::Velocity>) {
                const auto& vel = dataSource.GetVelocity();
                field_buffers[i].push_back(vel.x());
                field_buffers[i].push_back(vel.y());
                field_buffers[i].push_back(vel.z());
                // field_buffers[i].push_back(std::isinf(vel.x()) ? 0.0 : vel.x());
                // field_buffers[i].push_back(std::isinf(vel.y()) ? 0.0 : vel.y());
                // field_buffers[i].push_back(std::isinf(vel.z()) ? 0.0 : vel.z());
             } else if constexpr (std::is_same_v<std::decay_t<decltype(src_type)>, source::ShearStress>) {
                double val = dataSource.GetShearStress();
                field_buffers[i].push_back(val);
                // field_buffers[i].push_back(std::isinf(val) ? -1.0 : val);
             }
          }, outputSpec.fields[i].src);
        }
      }
    }
    
    // write datasets once gathered
    WriteDataset(group_id, "geometry", coords_buffer, global_site_count, local_site_count, 3, comms, mpi_comm);
    for (size_t i = 0; i < outputSpec.fields.size(); ++i) {
        const auto& field_spec = outputSpec.fields[i];
        WriteDataset(group_id, field_spec.name, field_buffers[i], global_site_count, local_site_count, GetFieldLength(field_spec.src), comms, mpi_comm);
    }

    H5_CHECK(H5Gclose(group_id));
  }
}
#endif