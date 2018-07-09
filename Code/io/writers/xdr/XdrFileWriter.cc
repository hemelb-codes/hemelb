
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cstdio>
#include "io/writers/xdr/XdrFileWriter.h"
#include "Exception.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

        // Implement a constructor that opens the file and creates the Xdr
        // object to write to it.
        XdrFileWriter::XdrFileWriter(const std::string& fileName, const std::string& mode)
        {
          myFile = std::fopen(fileName.c_str(), mode.c_str());
          if (myFile == NULL)
          {
            throw Exception() << "Failed to open file '" << fileName << "'";
          }
          xdrstdio_create(&mXdr, myFile, XDR_ENCODE);
        }

        // A destructor that ends the work of the Xdr object (including a
        // flush to file, so the order is important here), then frees the
        // memory and closes the file.
        XdrFileWriter::~XdrFileWriter()
        {
          xdr_destroy(&mXdr);
          std::fclose(myFile);
        }

      } // namespace xdr

    } // namespace writers

  }
}
