#include "io/hdf5/H5Error.h"
#include "log/Logger.h"

#include <hdf5.h>

namespace hemelb
{
  namespace io
  {

    namespace hdf5
    {

      H5Error::H5Error(const std::string & function, herr_t error,
                       const std::string & file, unsigned int line) :
          function(function), error(error), file(file), line(line) {
        operator<<("HDF5 error code ");
        operator<<(error);
        operator<<(" in ");
        operator<<(function);
        operator<<(" (");
        operator<<(file);
        operator<<(":");
        operator<<(line);
        operator<<(")");
      }

      void H5FileDeleter(hid_t * id) {
        herr_t error = H5Fclose(*id);
        delete id;
        if (error < 0)
        {
          // Log but don't throw
          hemelb::log::Logger::Log<hemelb::log::Error, hemelb::log::Singleton>(
              "HemeLB Error in %s (%s:%d): Failed to close HDF5 file (error %d)",
              __func__, __FILE__, __LINE__, error);
        }
      }

      void H5GroupDeleter(hid_t * id) {
        herr_t error = H5Gclose(*id);
        delete id;
        if (error < 0)
        {
          // Log but don't throw
          hemelb::log::Logger::Log<hemelb::log::Error, hemelb::log::Singleton>(
              "HemeLB Error in %s (%s:%d): Failed to close HDF5 group (error %d)",
              __func__, __FILE__, __LINE__, error);
        }
      }

      void H5DataSetDeleter(hid_t * id) {
        herr_t error = H5Dclose(*id);
        delete id;
        if (error < 0)
        {
          // Log but don't throw
          hemelb::log::Logger::Log<hemelb::log::Error, hemelb::log::Singleton>(
              "HemeLB Error in %s (%s:%d): Failed to close HDF5 dataset (error "
              "%d)", __func__, __FILE__, __LINE__, error);
        }
      }

      void H5DataSpaceDeleter(hid_t * id) {
        herr_t error = H5Sclose(*id);
        delete id;
        if (error < 0)
        {
          // Log but don't throw
          hemelb::log::Logger::Log<hemelb::log::Error, hemelb::log::Singleton>(
              "HemeLB Error in %s (%s:%d): Failed to close HDF5 data space "
              "(error %d)", __func__, __FILE__, __LINE__, error);
        }
      }

      void H5AttributeDeleter(hid_t * id) {
        herr_t error = H5Aclose(*id);
        delete id;
        if (error < 0)
        {
          // Log but don't throw
          hemelb::log::Logger::Log<hemelb::log::Error, hemelb::log::Singleton>(
              "HemeLB Error in %s (%s:%d): Failed to close HDF5 attribute "
              "(error %d)", __func__, __FILE__, __LINE__, error);
        }
      }

      void H5PropertyListDeleter(hid_t * id) {
        herr_t error = H5Pclose(*id);
        delete id;
        if (error < 0)
        {
          // Log but don't throw
          hemelb::log::Logger::Log<hemelb::log::Error, hemelb::log::Singleton>(
              "HemeLB Error in %s (%s:%d): Failed to close HDF5 property list "
              "(error %d)", __func__, __FILE__, __LINE__, error);
        }
      }
    }
  }
}
