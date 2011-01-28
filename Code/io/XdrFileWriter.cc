#include <stdio.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "io/XdrFileWriter.h"

namespace hemelb
{
  namespace io
  {

    // Implement a constructor that opens the file and creates the Xdr
    // object to write to it.
    XdrFileWriter::XdrFileWriter(std::string fileName)
    {
      myFile = fopen(fileName.c_str(), "w");
      mXdr = new XDR;
      xdrstdio_create(mXdr, myFile, XDR_ENCODE);
    }

    // A destructor that ends the work of the Xdr object (including a
    // flush to file, so the order is important here), then frees the
    // memory and closes the file.
    XdrFileWriter::~XdrFileWriter()
    {
      xdr_destroy(mXdr);
      delete mXdr;
      fclose(myFile);
    }

  }
}
