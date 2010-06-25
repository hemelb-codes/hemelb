#include "xdrFileWriter.h"

// Implement a constructor that opens the file and creates the Xdr object to write to it.
XdrFileWriter ::XdrFileWriter(char* fileName)
{
  myFile = fopen (fileName, "w");
  xdrstdio_create (&myXdr, myFile, XDR_ENCODE);
}

// A destructor that ends the work of the Xdr object (including a flush to file, so the order is
// important here), then closes the file.
XdrFileWriter ::~XdrFileWriter()
{
  xdr_destroy (&myXdr);
  fclose (myFile);
}
