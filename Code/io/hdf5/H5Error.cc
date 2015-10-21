#include "io/hdf5/H5Error.h"

#include <cstdlib>
#include <cstdio>
#include <hdf5.h>

#include "log/Logger.h"

static herr_t error_walker(unsigned int n, const H5E_error2_t * err_desc,
                           void * client_data)
{
  ssize_t class_size = H5Eget_class_name(err_desc->cls_id, NULL, 0);
  if (class_size < 0)
    return class_size;
  char class_name[class_size];
  if (class_size == 0)
    class_name[0] = 0;
  else if ((class_size = H5Eget_class_name(err_desc->cls_id, class_name,
                                           class_size)) < 0)
    return class_size;

  char ** message = (char **)client_data;
  int size = std::snprintf(*message, 0,
                           "Error (class %s, number %d.%d) in %s (%s:%u): %s\n",
                           class_name, err_desc->maj_num, err_desc->min_num,
                           err_desc->func_name, err_desc->file_name,
                           err_desc->line, err_desc->desc);
  if (size < 0)
    return size;
  if ((*message = (char *)std::realloc(message, (std::size_t)size + 1)) == NULL)
    return -1;
  if (std::snprintf(*message, (std::size_t)size,
                    "Error (class %s, number %d.%d) in %s (%s:%u): %s\n",
                    class_name, err_desc->maj_num, err_desc->min_num,
                    err_desc->func_name, err_desc->file_name,
                    err_desc->line, err_desc->desc) < 0)
    return -1;
  return 0;
}

namespace hemelb
{
  namespace io
  {

    namespace hdf5
    {

      H5Error::H5Error()
      {
        char * error_message = NULL;
        if (H5Ewalk2(H5E_DEFAULT, H5E_WALK_UPWARD, error_walker,
                     (void *)&error_message) < 0)
        {
          hemelb::log::Logger::Log<hemelb::log::Debug,
                                   hemelb::log::OnePerCore>("Closing HDF5 identifier failed: unknown type");
        }
        wat.assign(error_message);
        std::free(error_message);
      }
    }
  }
}
