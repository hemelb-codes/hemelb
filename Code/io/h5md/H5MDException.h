#ifndef HEMELB_IO_H5MD_H5MDEXCEPTION_H
#define HEMELB_IO_H5MD_H5MDEXCEPTION_H

#include "Exception.h"

namespace hemelb
{
  namespace io
  {
    namespace h5md
    {

      class H5MDException : public ::hemelb::Exception {
        public:
          H5MDException() {}
          H5MDException(const std::string & what) : reason(what) {}
          const std::string & what() { return reason; }
        private:
          std::string reason;
      };

    }
  }
}

#endif  // HEMELB_IO_H5MD_H5MDEXCEPTION_H
