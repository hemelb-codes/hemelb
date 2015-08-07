#include "io/h5md/TimeData.h"

#include <cstdlib>
#include <string>
#include <H5Cpp.h>

namespace hemelb
{
  namespace io
  {
    namespace h5md
    {
      TimeData TimeData::CreateExplicit(const H5::Group & location,
                                     const std::string & name,
                                     int rank, const hsize_t * dims,
                                     const hsize_t * max_dims = dims,
                                     const hsize_t * chunk_dims = dims,
                                     const H5::DataType & value_type = H5::PredType::NATIVE_DOUBLE,
                                     const H5::DataType & step_type = H5::PredType::NATIVE_ULONG,
                                     const H5::DataType & time_type = H5::PredType::NATIVE_DOUBLE,
                                     hsize_t timestep_chunk = 1,
                                     bool trackTime = false)
      {
        if (!value_type.detectClass(H5T_ENUM) &&
            !value_type.detectClass(H5T_INTEGER) &&
            !value_type.detectClass(H5T_FLOAT) &&
            !value_type.detectClass(H5T_STRING))
        {
          throw H5MDException("Data type must be one of Enumeration, Integer, Float or String");
        }
        if (!step_type.detectClass(H5T_INTEGER))
        {
          throw H5MDException("Time step type must be Integer");
        }
        if (!time_type.detectClass(H5T_FLOAT))
        {
          throw H5MDException("Time type must be Float");
        }

        // Create a group with the specified name in the specified location
        const H5::Group group = location.createGroup(name);

        // Expand the dimensions arrays to include the unlimited time dimension
        hsize_t expanded_dims[rank + 1];
        expanded_dims[0] = H5S_UNLIMITED;
        std::copy_n(dims, rank, &expanded_dims[1]);
        hsize_t expanded_max_dims[rank + 1];
        expanded_max_dims[0] = H5S_UNLIMITED;
        std::copy_n(max_dims, rank, &expanded_max_dims[1]);
        hsize_t expanded_chunk_dims[rank + 1];
        expanded_chunk_dims[0] = timestep_chunk;
        std::copy_n(chunk_dims, rank, &expanded_chunk_dims[1]);
        rank = rank + 1;

        // Create the value dataset
        const H5::DataSpace value_space(rank, expanded_dims, expanded_max_dims);
        const H5::DSetCreatPropList create_plist;
        create_plist.setLayout(H5D_CHUNKED);
        create_plist.setChunk(rank, chunk_dims);
        const H5::DataSet value = group.createDataSet("value",
                                                      value_type,
                                                      value_space,
                                                      create_plist);

        // Create the timestep datasets
        const hsize_t timestep_dims[] = { H5S_UNLIMITED };
        const H5::DataSpace timestep_space(1, timestep_dims);
        create_plist.setChunk(1, &timestep_chunk);
        const H5::DataSet step = group.createDataSet("step",
                                                     step_type,
                                                     timestep_space,
                                                     create_plist);
        H5::DataSet time;
        if (trackTime)
        {
          time = group.createDataSet("time", time_type, timestep_space,
                                     create_plist);
        }

        return TimeData(group, name, value, step, time);
      }

      TimeData TimeData::CreateFixed(const H5::Group & location,
                                  const std::string & name,
                                  int rank, const hsize_t * dims,
                                  const hsize_t * max_dims = dims,
                                  const hsize_t * chunk_dims = dims,
                                  unsigned long step_offset = 0, double time_offset = 0.0,
                                  const H5::DataType & value_type = H5::PredType::NATIVE_DOUBLE,
                                  const H5::DataType & step_type = H5::PredType::NATIVE_ULONG,
                                  const H5::DataType & time_type = H5::PredType::NATIVE_DOUBLE,
                                  hsize_t timestep_chunk = 1,
                                  bool trackTime = false)
      {
        if (!value_type.detectClass(H5T_ENUM) &&
            !value_type.detectClass(H5T_INTEGER) &&
            !value_type.detectClass(H5T_FLOAT) &&
            !value_type.detectClass(H5T_STRING))
        {
          throw H5MDException("Data type must be one of Enumeration, Integer, Float or String");
        }
        if (!step_type.detectClass(H5T_INTEGER))
        {
          throw H5MDException("Time step type must be Integer");
        }
        if (!time_type.detectClass(H5T_FLOAT))
        {
          throw H5MDException("Time type must be Float");
        }

        // Create a group with the specified name in the specified location
        const H5::Group group = location.createGroup(name);

        // Expand the dimensions arrays to include the unlimited time dimension
        hsize_t expanded_dims[rank + 1];
        expanded_dims[0] = H5S_UNLIMITED;
        std::copy_n(dims, rank, &expanded_dims[1]);
        hsize_t expanded_max_dims[rank + 1];
        expanded_max_dims[0] = H5S_UNLIMITED;
        std::copy_n(max_dims, rank, &expanded_max_dims[1]);
        hsize_t expanded_chunk_dims[rank + 1];
        expanded_chunk_dims[0] = timestep_chunk;
        std::copy_n(chunk_dims, rank, &expanded_chunk_dims[1]);
        rank = rank + 1;

        // Create the value dataset
        const H5::DataSpace value_space(rank, expanded_dims, expanded_max_dims);
        const H5::DSetCreatPropList create_plist;
        create_plist.setLayout(H5D_CHUNKED);
        create_plist.setChunk(rank, chunk_dims);
        const H5::DataSet value = group.createDataSet("value",
                                                      value_type,
                                                      value_space,
                                                      create_plist);

        // Create the timestep datasets
        const H5::DataSpace timestep_space;
        const H5::DataSet step = group.createDataSet("step", step_type, timestep_space);
        if (step_offset > 0)
        {
          H5::Attribute attr = step.createAttribute("offset", step_type, timestep_space);
          attr.write(step_type, &step_offset);
        }
        H5::DataSet time;
        if (trackTime)
        {
          time = group.createDataSet("time", time_type, timestep_space);
          if (time_offset > 0.0)
          {
            H5::Attribute attr = time.createAttribute("offset", time_type, timestep_space);
            attr.write(time_type, &time_offset);
          }
        }

        return TimeData(group, name, value, step, time);
      }

      TimeData TimeData::Open(const H5::Group & location,
                           const std::string & name) {
        // Open the group
        const H5::Group group = location.openGroup(name);

        // Open the value, step and time datasets
        const H5::DataSet value = group.openDataSet("value");
        const H5::DataSet step = group.openDataSet("step");
        const H5::DataSet time = group.openDataSet("time");

        return TimeData(group, name, value, step, time);
      }

      bool TimeData::IsStepExplicit() const
      {
        return step.getSpace().getSimpleExtentType() != H5S_SCALAR;
      }

      bool TimeData::HasTime() const
      {
        try {
          group.openDataSet("time");
        } catch (H5::Exception &) {
          return false;
        }
        return true;
      }

    }
  }
}
