#ifndef HEMELB_IO_H5MD_H5MDERROR_H
#define HEMELB_IO_H5MD_H5MDERROR_H

#include <exception>
#include <string>

namespace hemelb
{
  namespace io
  {
    namespace h5md
    {

      class H5MDError : public std::exception
      {
        public:
          H5MDError()
          {
          }
          H5MDError(const std::string & what) :
              wat(what)
          {
          }
          virtual const char * what() const noexcept
          {
            return wat.c_str();
          }
        private:
          std::string wat;
      };

    }
  }
}

#endif  // HEMELB_IO_H5MD_H5MDERROR_H
